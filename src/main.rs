use core::{panic};
use std::{io::{BufRead, Write, BufReader}, fs::File, string, collections::HashMap};

use bio::io::fasta::Record;
use clap::{Parser, Arg};
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IndexedParallelIterator};

mod reference;
mod find_cdr3;
mod find_v;

#[derive(Parser,Debug)]
struct Args {
    #[arg(short,long)]
    input_fasta: String,
    #[arg(short,long)]
    reference_fasta: String,
    #[arg(short,long)]
    output_csv: String,

    #[arg(short,long, default_value_t = 10000)]
    sample_size: usize,
    #[arg(short,long, default_value_t = 500000)]
    parallel_chunk_size: usize,
    #[arg(short,long)]
    nondeterministic: bool,
    #[arg(short,long, default_value_t = 10)]
    edit_dist: u8,
}

mod output {
    use std::fs::File;
    use std::io::Write;
    use crate::find_cdr3::OutputRecord;

    pub fn print_header(file: &mut File) {
        write!(file, "id,sequence,cdr3_sequence\n")
        .expect("Couldn't write header line to output!");
    }
    
    pub fn print_one(output_file: &mut File, output_record: OutputRecord) {
        write!(output_file, "{}\n", output_record)
            .expect("Couldn't write line to output!");
    }
}

/// The main function with all the stuff in it
fn main() {
    // get the arguments from the command line
    let args = Args::parse();

    // collect all the reference seqs
    let reference_seqs = reference::parse_reference(&args.reference_fasta);
    
    // open the files to read and write from
    let records = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(&args.input_fasta)
                .expect("Couldn't open input!")))
                .records();
    let mut output_csv_file = File::create(args.output_csv)
        .expect("Couldn't create output file!");
    output::print_header(&mut output_csv_file);

    // optimise the reference seqs list
    let optimised_reference_seqs = find_v::optimise_refs(
        &reference_seqs, 
        bio::io::fasta::Reader::new(
            BufReader::new(
                File::open(&args.input_fasta)
                    .expect("Couldn't open input!")))
            .records(), args.sample_size, args.parallel_chunk_size, args.edit_dist);

    // go through each result and print to output while you go 
    // to avoid collecting data in memory
    
    for bunch in &records.chunks(args.parallel_chunk_size) {
        let mut outs = Vec::new();
        
        bunch.collect_vec().into_par_iter().map(|result| {
            match result {
                Err(_) => panic!("Bad record!"),
                Ok(record) => if !args.nondeterministic {
                    find_cdr3::det::parse_one_input_par(
                        record, &optimised_reference_seqs, args.edit_dist)
                } else {
                    find_cdr3::nondet::parse_one_input_par(
                        record, &optimised_reference_seqs, args.edit_dist)
                }

            }
        }).collect_into_vec(&mut outs);

        for output in outs {
            output::print_one(&mut output_csv_file, output)
        }
    }

    // for result in records {
    //     match result {
    //         Err(_) => panic!("Bad record!"),
    //         Ok(record) => {
    //             output::print_one(
    //                 &mut output_csv_file,
    //                 find_cdr3::det::parse_one_input(record, &optimised_reference_seqs, 20)
    //             );
    //         }
    //     }
    // }
}

fn mot_main() {
    // get the arguments from the command line
    let args = Args::parse();

    // collect all the reference seqs
    let reference_seqs = reference::parse_reference(&args.reference_fasta);
    
    // open the files to read and write from
    let records = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(args.input_fasta)
                .expect("Couldn't open input!")))
                .records();
    let mut output_csv_file = File::create(args.output_csv)
        .expect("Couldn't create output file!");
    output::print_header(&mut output_csv_file);

}

struct Feature<'a> {
    seq: &'a str,
    
    start: usize,
    end: usize,
}

struct ImportantExcerpt<'a> {
    part: &'a str
}

/// Quick function just for fun
/// # Examples
/// ```
/// let s = String::from("hello");
/// let b = quick(s);
/// ```
fn quick<'a>(inp: &'a String) -> ImportantExcerpt<'a> {
    let impo = ImportantExcerpt {
        part: &inp[1..]
    };

    impo
}