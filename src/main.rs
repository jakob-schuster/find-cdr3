use core::{panic};
use std::{io::{BufRead, Write}, fs::File, string};

use bio::{pattern_matching::myers::{MyersBuilder}, io::fasta::Record, alphabets::protein, alignment::{}};
use clap::{Parser, Arg};

mod reference;
mod find_cdr3;

#[derive(Parser,Debug)]
struct Args {
    #[arg(short,long)]
    input_fasta: String,
    #[arg(short,long)]
    reference_fasta: String,
    #[arg(short,long)]
    output_csv: String,
}



fn print_output(output_csv: String, output_records: &Vec<find_cdr3::OutputRecord>) {
    let mut output_csv_file = File::create(output_csv)
        .expect("Couldn't create output file!");

    write!(output_csv_file, "id,sequence,cdr3_sequence\n")
        .expect("Couldn't write header line to output!");
    for record in output_records {
        write!(output_csv_file, "{}\n", record)
            .expect("Couldn't write line to output!");
    }
}

fn main() {
    let args = Args::parse();

    let reference_seqs = reference::parse_reference(args.reference_fasta);
    let output_records = find_cdr3::parse_input(args.input_fasta, &reference_seqs, 8);

    print_output(args.output_csv, &output_records);
}