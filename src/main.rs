use std::any::{Any, TypeId};

use clap::Parser;
use itertools::Itertools;

use crate::{reference::EfficientGuides, find::{StructureResult, RefResult}};

mod find;
mod reference;

#[derive(Parser,Debug)]
pub struct Args {
    #[arg(short, long, default_value_t = String::from("stdin"))]
    input_fastq: String,

    #[arg(short,long)]
    reference_tsv: String,

    #[arg(short,long, default_value_t = String::from("stdout"))]
    output: String,

    #[arg(short, long, default_value_t = String::from("GTTCACTGCCGTATAGGCAG"))]
    cys4: String,

    #[arg(short, long, default_value_t = String::from("GAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"))]
    scaffold: String,

    /// Edit distance used for reference sequences.
    #[arg(short,long, default_value_t = 10)]
    edit_dist: u8,
}

fn _seq_to_string(seq: &[u8]) -> String {
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

fn main() {
    println!("Started..");

    let arg = Args::parse();
    let reference = reference::Ref::new(&arg);

    println!("Parsed reference..");

    let records = bio::io::fastq::Reader::new(input::reader(&arg))
        .records();
    
    let mut out = vec![];

    let efficient_guides = EfficientGuides::new(&reference.guides);
    println!("Produced efficient reference..");

    for result in records {
        match result {
            Ok(record) => {
                out.push(find::structure_classify(record.seq(), &reference, &efficient_guides));
            },
            Err(_) => panic!("Bad record!"),
        }
    }

    print_stats(&out);
}

fn print_stats(out: &Vec<StructureResult>) {
    let total = out.clone().len();
    
    let well_structured = out.clone().into_iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(_))
    }).collect_vec().len();

    let valid = out.clone().into_iter().filter(|r| {
        matches!(r, StructureResult::WellStructured(RefResult::Valid(_)))
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