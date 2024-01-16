use core::{panic};
use std::{io::{BufRead, Write, BufReader, stdout, stdin, Stdin, BufWriter}, fs::File, string, collections::HashMap, path::Path, ffi::OsStr};

use bio::{io::fasta::Record, pattern_matching::myers::Myers, alignment::{Alignment, AlignmentOperation}};
use clap::{Parser, Arg, ValueHint};
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IndexedParallelIterator};

mod reference;
mod find_cdr3;
mod find_v;

#[derive(Parser,Debug)]
pub struct Args {
    /// Input fasta file of immunoglobin sequences. Reads from stdin if none is provided.
    #[arg(short, long, default_value_t = String::from("stdin"))]
    input_fasta: String,
    /// Reference fasta of different V-gene sequences. The Cys codon near the end of the V-gene sequence marks the start of the CDR3.
    #[arg(short,long)]
    reference_fasta: String,
    /// Output TSV file of reads and their CDR3s. Writes to stdout if none is provided.
    #[arg(short,long, default_value_t = String::from("stdout"))]
    output_tsv: String,

    /// Size of pre-processing sample, pre-processed to find the most common reference sequences. Not used currently.
    #[arg(long, default_value_t = 10000)]
    sample_size: usize,
    /// Number of top V-gene sequences to actually search for. The less you use, the quicker the program, but also the less accurate.
    #[arg(long, default_value_t = 20)]
    reference_size: usize,
    
    /// Size of each chunk of reads to process in parallel.
    #[arg(long, short='c',default_value_t = 500000)]
    parallel_chunk_size: usize,
    /// Experimental. Use the non-deterministic approach to searching. No performance gains yet.
    #[arg(short,long,hide=true)]
    nondeterministic: bool,
    /// Edit distance used for reference sequences.
    #[arg(short,long, default_value_t = 20)]
    edit_dist: u8,
    /// FR4 region. Located just after the CDR3.
    #[arg(short,long, default_value_t = String::from("TGGGGCAAAGGGACCCAGGTCAC"))]
    fr4: String,
    /// Omit the header from the output. Useful if you're piping into another program to do more processing.
    #[arg(short='l', long,default_value_t = false)]
    headerless: bool
}

mod input {
    use std::{io::{BufRead, BufReader, stdin}, path::Path, fs::File, ffi::OsStr};

    use crate::Args;

    /// Checks the arguments, and either opens a file or reads from stdin
    pub fn reader(args: &Args) -> Box<dyn BufRead> {
        if args.input_fasta.eq("stdin") {
            // just read straight from stdin
            Box::new(BufReader::new(stdin()))
        } else {
            let path = Path::new(&args.input_fasta);
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
    use crate::find_cdr3::OutputRecord;

    pub fn print_header<T: Write>(output: &mut T) {
        write!(output, "id\tsequence\tcdr3_sequence\n")
        .expect("Couldn't write header line to output!");
    }
    
    pub fn print_one<T: Write>(output: &mut T, output_record: OutputRecord) {
        write!(output, "{}\n", output_record)
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

/// The main function with all the stuff in it
fn main() {
    // get the arguments from the command line
    let args = Args::parse();

    // collect all the reference seqs
    let reference_seqs = reference::parse_reference(&args.reference_fasta);

    // open the files to read and write from
    let records = bio::io::fasta::Reader::new(input::reader(&args))
        .records();
    
    // optimise the reference seqs list
    let (optimised_reference_seqs, _rest) = find_v::optimise_refs(
        &reference_seqs, 
        bio::io::fasta::Reader::new(input::reader(&args)).records(),
        args.sample_size, args.parallel_chunk_size, args.edit_dist, args.reference_size);
    // let optimised_reference_seqs = reference_seqs;
    
    let mut output_csv_file = output::writer(&args);
    // print the header unless the user specifies not to
    if !args.headerless {
        output::print_header(&mut output_csv_file);
    }

    let fr4 = Myers::<u64>::new(args.fr4.as_bytes());

    // select a parsing function based on the nondeterminism parameter
    let parse_one = |record| if !args.nondeterministic {
        find_cdr3::det::parse_one_input_par(
            record, 
            &optimised_reference_seqs, 
            args.edit_dist, 
            args.fr4.as_bytes()
        )
    } else {
        find_cdr3::nondet::parse_one_input(
            record, 
            &optimised_reference_seqs, 
            args.edit_dist, 
            &fr4
        )
    };

    // go through each result and print to output while you go 
    // to avoid collecting data in memory
    for bunch in &records.chunks(args.parallel_chunk_size) {
        let mut outs = Vec::new();
        
        bunch.collect_vec().into_par_iter().map(|result| {
            match result {
                Err(_) => panic!("Bad record!"),
                Ok(record) => parse_one(record)
            }
        }).collect_into_vec(&mut outs);

        for output in outs {
            output::print_one(&mut output_csv_file, output)
        }
    }
}

fn to_char(op: &AlignmentOperation) -> char {
    match op {
        AlignmentOperation::Match => 'm',
        AlignmentOperation::Subst => 's',
        AlignmentOperation::Del => 'd',
        AlignmentOperation::Ins => 'i',
        AlignmentOperation::Xclip(_) => 'x',
        AlignmentOperation::Yclip(_) => 'y',
    }
}

fn find_index(ops: &[AlignmentOperation], ref_index: usize) -> Option<usize> {
    let mut seq_index = None;
    let (mut seq_i, mut ref_i) = (0, 0);

    for op in ops {
        if ref_i == ref_index {
            seq_index = Some(seq_i);
            // break;
        }

        match op {
            AlignmentOperation::Del => { seq_i += 1 },
            AlignmentOperation::Ins => { ref_i += 1 },
            AlignmentOperation::Match => { seq_i += 1; ref_i += 1 },
            AlignmentOperation::Subst => { seq_i += 1; ref_i += 1 },
            AlignmentOperation::Xclip(_) => panic!("Xclip"),
            AlignmentOperation::Yclip(_) => panic!("Yclip")
        }
    }
    
    seq_index
}

fn not_main() {
    let seq = b"
AAAATGGGAAGAGTCGT";

    let pattern = b"AGGTGGAAG";
    let tgg_end_pos = 6;

    let mut myers = Myers::<u64>::new(pattern);

    let mut matches = myers.find_all_lazy(seq, 4);
    let best_match = 
        matches.by_ref().min_by_key(|&(_, dist)| dist);

    match best_match {
        Some((end, dist)) => {
            let mut aln = Alignment::default();
            matches.alignment_at(end, &mut aln);
            let ops = aln.operations;
            
            let ops_string: String = ops.clone().into_iter().map(|op| to_char(&op)).collect();
            let index = find_index(&ops[..], tgg_end_pos);

            match index {
                Some(i) => {
                    println!("bit was {} in orig. now it's {}", 
                        String::from_utf8(pattern[tgg_end_pos..].to_vec()).expect("aa"), 
                        String::from_utf8(seq[i+aln.ystart..].to_vec()).expect("aa")
                    );
                    println!("ystart {} in {}", aln.ystart, ops_string);
                }
            None => todo!(),
            }


        },
        None => todo!(),
    }
    
    println!("up to here");
}