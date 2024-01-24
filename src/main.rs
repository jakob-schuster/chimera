use std::any::{Any, TypeId};

use bio::stats;
use clap::{Parser, builder::Str};
use itertools::Itertools;
use output::writer;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IndexedParallelIterator};

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
    #[arg(short, long, default_value_t = String::from("stdin"))]
    input_fastq: String,

    #[arg(short,long)]
    reference_tsv: String,

    #[arg(short,long, default_value_t = String::from("stdout"))]
    output_tsv: String,

    #[arg(short, long, default_value_t = String::from("GTTCACTGCCGTATAGGCAG"))]
    cys4: String,
                                           
    #[arg(short, long, default_value_t = String::from("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"))]
    scaffold: String,

    /// Edit distance used for reference sequences.
    #[arg(short,long, default_value_t = 0.25)]
    error_rate: f32,
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
        write!(output, "id\tsequence\tcdr3_sequence\n")
        .expect("Couldn't write header line to output!");
    }
    
    pub fn print_one<T: Write>(output: &mut T, output_record: (String, StructureResult)) {
        write!(output, "{}\t{:?}\n", output_record.0, output_record.1)
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

fn main() {
    // println!("Started..");

    let args = Args::parse();
    let reference = reference::Ref::new(&args);

    // println!("Parsed reference..");

    let records = bio::io::fastq::Reader::new(input::reader(&args))
        .records();
    
    let mut writer = writer(&args);

    let efficient_guides = EfficientGuides::new(&reference.guides);
    // println!("Produced efficient reference..");

    let mut out = Vec::new();
    // for result in records {
    //     match result {
    //         Ok(record) => {out.push((String::from(record.id()), find::structure_classify(record.seq(), &reference, &efficient_guides)))},
    //         Err(_) => panic!("Bad record!"),
    //     }
    // }

    for chunk in &records.chunks(100000) {
        let mut temp = Vec::new();
        chunk.collect_vec().into_par_iter().map(|result| -> (String, StructureResult) {match result {
            Ok(record) => (String::from(record.id()), find::structure_classify(record.seq(), &reference, &efficient_guides, args.error_rate)),
            Err(_) => panic!("Bad record!"),
        }}).collect_into_vec(&mut temp);

        out.append(&mut temp);
    }

    for (id, structure) in out.clone() {
        writer.write(format!("{}\t{:?}\n", id, structure).as_bytes());
    }

    let stats_out = out.into_iter().map(|(_, s)| s).collect();
    print_stats(&stats_out);

}

fn print_stats(out: &Vec<StructureResult>) {
    let total = out.clone().len();
    
    let well_structured = out.clone().into_iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(_))
    }).collect_vec().len();

    let valid_reads_map = out.clone().into_iter().counts();
    let zipped = valid_reads_map.keys().into_iter().zip(valid_reads_map.values()).sorted_by_key(|(a, b)| **b);
    for (k, v) in zipped {
        println!("{:?}\t{}", k, v);
    }

    let valid = out.clone().into_iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(RefResult::Valid(_, _)))
    }).collect_vec().len();

    let chimeric = out.clone().into_iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(RefResult::Chimera))
    }).collect_vec().len();

    let ambiguous = out.clone().into_iter().filter(|r| {
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