use std::{any::{Any, TypeId}, fs::File, io::BufWriter};

use bio::stats;
use clap::{Parser, builder::Str};
use itertools::Itertools;
use output::writer;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IndexedParallelIterator};
use reference::{FinalGuides, Ref};
use seq_io::{fastq::RefRecord, parallel::parallel_fastq};

use crate::{reference::EfficientGuides, find::{StructureResult, RefResult}};

mod find;
mod reference;

// #[derive(Debug, Parser)]
// struct Cli {
    
// }

// enum Commands {
    
// }

#[derive(Parser,Debug)]
pub struct Args {
    #[arg(short, long)]
    input_fastq: String,

    #[arg(short,long)]
    reference_tsv: String,

    #[arg(long, default_value_t = false)]
    print_stats: bool,
    
    #[arg(short, long)]
    output_chimera_fastq: String,

    #[arg(short, long)]
    output_valid_fastq: String,

    #[arg(short, long)]
    output_tsv: String,

    #[arg(short, long, default_value_t = String::from("GTTCACTGCCGTATAGGCAG"))]
    cys4: String,
                                           
    #[arg(short, long, default_value_t = String::from("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"))]
    // #[arg(short, long, default_value_t = String::from("GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"))]
    scaffold: String,

    /// Edit distance used for reference sequences.
    #[arg(short,long, default_value_t = 0.25)]
    error_rate: f32,

    #[arg(long, default_value_t = false)]
    careful: bool,

}

pub fn _seq_to_string(seq: &[u8]) -> String {
    String::from_utf8(seq.to_vec()).expect("Bad seq!")
}

mod input {
    use std::{io::{BufRead, BufReader, stdin}, path::Path, fs::File, ffi::OsStr};

    use crate::Args;

    /// Checks the arguments, and either opens a file or reads from stdin
    pub fn reader(args: &Args) -> Box<dyn BufRead> {
        if args.input_fastq.eq("stdin") {
            // just read straight from stdin
            Box::new(BufReader::new(stdin()))
        } else {
            let path = Path::new(&args.input_fastq);
            let file = match File::open(path) {
                Ok(file) => file,
                Err(_) => panic!("Couldn't open {}!", path.display()),
            };

            if path.extension() == Some(OsStr::new("gz")) {
                Box::new(BufReader::new(
                    flate2::read::MultiGzDecoder::new(file)))
            } else {
                Box::new(BufReader::new(file))
            }
        }
    }
}

mod output {
    use std::ffi::OsStr;
    use std::fs::File;
    use std::io::{Write, BufWriter, stdout};
    use std::path::Path;
    use crate::Args;
    use crate::find::StructureResult;

    pub fn print_header<T: Write>(output: &mut T) {
        writeln!(output, "id\tsequence\tcdr3_sequence")
        .expect("Couldn't write header line to output!");
    }
    
    pub fn print_one<T: Write>(output: &mut T, output_record: (String, StructureResult)) {
        writeln!(output, "{}\t{:?}", output_record.0, output_record.1)
            .expect("Couldn't write line to output!");
    }

    /// Checks the arguments, and either opens a file or writes to stdout
    pub fn writer(args: &Args) -> Box<dyn Write> {
        if args.output_tsv.eq("stdout") {
            // just read straight from stdin
            Box::new(BufWriter::new(stdout()))
        } else {
            let path = Path::new(&args.output_tsv);
            let file = match File::create(path) {
                Ok(file) => file,
                Err(_) => panic!("Couldn't open {}!", path.display()),
            };

            if path.extension() == Some(OsStr::new("gz")) {
                Box::new(BufWriter::new(
                    flate2::write::GzEncoder::new(file, flate2::Compression::default())))
            } else {
                Box::new(BufWriter::new(file))
            }
        }
    }
}

// fn main() {
//     // println!("Started..");

//     let args = Args::parse();
//     let reference = reference::Ref::new(&args);

//     // println!("Parsed reference..");

//     let records = bio::io::fastq::Reader::new(input::reader(&args))
//         .records();
    
//     let mut writer = writer(&args);

//     let efficient_guides = EfficientGuides::new(&reference.guides);
//     let final_guides = FinalGuides::new(&reference.guides);
//     // println!("Produced efficient reference..");

//     let mut out = Vec::new();

//     let classify = |seq: &[u8]| -> StructureResult {
//         if args.careful {
//             find::structure_classify_carefully(seq, &reference, &efficient_guides, args.error_rate)
//         } else {
//             find::structure_classify_quickly(seq, &reference, &final_guides, args.error_rate)
//         }
//     };

//     for chunk in &records.chunks(100000) {
//         let mut temp = Vec::new();
//         chunk.collect_vec().into_par_iter()
//             .map(|result| -> (String, StructureResult) {
//                 match result {
//                     Ok(record) => (String::from(record.id()), classify(record.seq())),
//                     Err(_) => panic!("Bad record!"),
//                 }
//             }
//         ).collect_into_vec(&mut temp);

//         out.append(&mut temp);
//     }

//     for (id, structure) in out.clone() {
//         let _ = writer.write(format!("{}\t{:?}\n", id, structure).as_bytes());
//     }

//     let stats_out = out.into_iter().map(|(_, s)| s).collect_vec();
//     print_stats(&stats_out);

// }


struct OutStats {
    total: u32,
    well_structured: u32,
    chimeric: u32,
    valid: u32,
    ambiguous: u32
}

impl OutStats {
    fn new() -> Self {
        OutStats {
            total: 0, 
            well_structured: 0, 
            chimeric: 0, 
            valid: 0, 
            ambiguous: 0 
        }
    }

    fn add(&mut self, result: &StructureResult) {
        self.total += 1;
        match result {
            StructureResult::WellStructured(w) => {
                self.well_structured += 1;
                match w {
                    RefResult::Valid(_, _) => { self.valid += 1 },
                    RefResult::Chimera => { self.chimeric += 1 },
                    RefResult::Ambiguous => { self.ambiguous += 1 },
                }
            },
            StructureResult::BadlyStructured => { },
        }
    }

    fn print_stats(&self) {
        println!("well-structured (scaffold - cys4 - scaffold): {} / {} = {}%", 
            self.well_structured, self.total, (self.well_structured as f32) / (self.total as f32) * 100.0);
        println!("valid (spacer, extension, nicking): {} / {} = {}% ({}% of well-structured reads)", 
            self.valid, self.total, (self.valid as f32) / (self.total as f32) * 100.0, (self.valid as f32) / (self.well_structured as f32) * 100.0);
        println!("chimeric (spacer, extension, nicking): {} / {} = {}% ({}% of well-structured reads)", 
                self.chimeric, self.total, (self.chimeric as f32) / (self.total as f32) * 100.0, (self.chimeric as f32) / (self.well_structured as f32) * 100.0);
        println!("ambiguous (spacer, extension, nicking): {} / {} = {}% ({}% of well-structured reads)", 
            self.ambiguous, self.total, (self.ambiguous as f32) / (self.total as f32) * 100.0, (self.ambiguous as f32) / (self.well_structured as f32) * 100.0);
    }
}

fn main() {
    use seq_io::fastq::{Reader,Record};

    let args = Args::parse();
    let reference = reference::Ref::new(&args);

    let reader = Reader::from_path(&args.input_fastq).unwrap();

    let mut writer = writer(&args);
    let mut valid_fastq = BufWriter::new(File::create(args.output_valid_fastq).unwrap());
    let mut chimera_fastq = BufWriter::new(File::create(args.output_chimera_fastq).unwrap());

    let efficient_guides = EfficientGuides::new(&reference.guides);
    let final_guides = FinalGuides::new(&reference.guides);


    let classify = |seq: &[u8]| -> StructureResult {
        if args.careful {
            find::structure_classify_carefully(seq, &reference, &efficient_guides, args.error_rate)
        } else {
            find::structure_classify_quickly(seq, &reference, &final_guides, args.error_rate)
        }
    };

    let f = |record: RefRecord<'_>, out: &mut StructureResult| {

    };

    // let mut out_vec = Vec::new();

    let mut out_stats = OutStats::new();
    let _ = parallel_fastq(reader, 32, 10000, |record, out| {
        *out = classify(record.seq());
    }, |record, out| {
        if let StructureResult::WellStructured(w) = out {
            match w {
                // write down the valid ones
                RefResult::Valid(_, _) => {
                    record.write(&mut valid_fastq).unwrap(); 
                },
                // write down the chimeras
                RefResult::Chimera => { 
                    record.write(&mut chimera_fastq).unwrap(); 
                },
                _ => {}
            }

            output::print_one(&mut writer, (String::from(record.id().unwrap()), out.to_owned()));
        }

        // keep them all in a big vec
        out_stats.add(out);
        // out_vec.push(out.clone());

        None::<()>
    });

    out_stats.print_stats()
    // print_stats(&out_vec);
}

fn print_stats(out: &[StructureResult]) {
    let total = out.len();
    
    let well_structured = out.iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(_))
    }).collect_vec().len();

    let valid_reads_map = out.iter().counts();
    let zipped = valid_reads_map.keys().zip(valid_reads_map.values()).sorted_by_key(|(_, b)| **b);
    for (k, v) in zipped {
        println!("{:?}\t{}", k, v);
    }

    let valid = out.iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(RefResult::Valid(_, _)))
    }).collect_vec().len();

    let chimeric = out.iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(RefResult::Chimera))
    }).collect_vec().len();

    let ambiguous = out.iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(RefResult::Ambiguous))
    }).collect_vec().len();

    println!("well-structured (scaffold - cys4 - scaffold): {} / {} = {}%", 
        well_structured, total, (well_structured as f32) / (total as f32) * 100.0);
    println!("valid (spacer, extension, nicking): {} / {} = {}% ({}% of well-structured reads)", 
        valid, total, (valid as f32) / (total as f32) * 100.0, (valid as f32) / (well_structured as f32) * 100.0);
    println!("chimeric (spacer, extension, nicking): {} / {} = {}% ({}% of well-structured reads)", 
            chimeric, total, (chimeric as f32) / (total as f32) * 100.0, (chimeric as f32) / (well_structured as f32) * 100.0);
    println!("ambiguous (spacer, extension, nicking): {} / {} = {}% ({}% of well-structured reads)", 
        ambiguous, total, (ambiguous as f32) / (total as f32) * 100.0, (ambiguous as f32) / (well_structured as f32) * 100.0);
}