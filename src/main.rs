use clap::Parser;

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
    let arg = Args::parse();
    let reference = reference::Ref::new(&arg);

    // let seq = b"GGCTTGTGGAAAGGACGAAAACCGCTGTCTCTTCCCCACAGAGGGTTTAAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCGACTTGAAAAAGTCGCACCGAGTCGGTGCGCTTCAAGATTTACTTTTGGCAATCGTCCTCTGTGGGGAAGTGCGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAAGTTGACTGCCGTATAGGCTGGCTCCTTCAAGAATTAGTTTGTTTTAGAGCTTGAAAATGAAAGTTTACATAGGGGTTGTCCGTTTTCAATTTTAAAACGTG";
    // let seq = 
    //     b"GTACTTGTGGAAAGGACGAAACACCGACTTTGAGTACTTGGAGATCGTTTAAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCGACTTGAAAAAGTCGCACCGAGTCGGTGCTCTCCAGTTGTCGGATCTCCAAGTACTCACGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAAGTTCACTGCCGTATAGGCAGCCGCTTGTGTCTCCAGTTGTGTTATAGTGCTAGAAATAGCAAGTTACAATAAGGCTAGTCGCTTTTCAACTTGAAAAAGTGCCACCGACTGGG";

    let records = bio::io::fastq::Reader::new(input::reader(&arg))
        .records();
    
    let mut out = vec![];
    let (mut ok, mut good, mut chimeric, mut total) = (0, 0, 0, 0);
    for result in records {
        match result {
            Ok(record) => {

                if let Some((spacer_seq, extension_seq, nicking_seq)) = 
                    find::break_into_regions(record.seq(), &reference) {
                    ok += 1;

                    if find::chimeric(
                        spacer_seq, 
                        extension_seq, 
                        nicking_seq, 
                        &reference
                    ) { 
                        chimeric += 1 
                    }
                }

                if let Some(name) = find::find(record.seq(), &reference) {
                    out.push((String::from(record.id()), name));

                    good += 1;
                }

                total += 1;
            },
            Err(_) => panic!("Bad record!"),
        }
    }

    println!("found scaffold and cys4 structure: {} / {} = {}", 
        ok, total, (ok as f32) / (total as f32) * 100.0);
    println!("found triple of (spacer, extension, nicking): {} / {} = {} ({}% of well-structured reads)", 
        good, total, (good as f32) / (total as f32) * 100.0, (good as f32) / (ok as f32) * 100.0);
    println!("chimeric (spacer, extension, nicking): {} / {} = {} ({}% of well-structured reads)", 
        chimeric, total, (chimeric as f32) / (total as f32) * 100.0, (chimeric as f32) / (ok as f32) * 100.0);
}