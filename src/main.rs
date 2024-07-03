use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use bio::{pattern_matching::myers::Myers, alignment::AlignmentOperation};
use clap::builder::Str;
use clap::Parser;
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use reference::RefV;
use seq_io::parallel::parallel_fasta;

mod reference;
mod find_cdr3;
mod find_v;

#[derive(Parser,Debug)]
pub struct Args {
    /// Input fasta file of immunoglobin sequences. Reads from stdin if none is provided.
    #[arg(short, long, default_value_t = String::from("stdin"))]
    input_reads: String,

    /// Reference fasta of different V-gene sequences. The Cys codon near the end of the V-gene sequence marks the start of the CDR3.
    #[arg(short, long)]
    reference_fasta: String,
    /// Output TSV file of reads and their CDR3s. Writes to stdout if none is provided.
    #[arg(short, long, default_value_t = String::from("stdout"))]
    output_tsv: String,
    /// Size of pre-processing sample, pre-processed to find the most common reference sequences.
    #[arg(long, default_value_t = 10000)]
    sample_size: usize,
    
    /// Number of top V-gene sequences to actually search for. The less you use, the quicker the program, but also the less accurate.
    #[arg(long, default_value_t = 5)]
    reference_size: usize,
    
    /// Size of each chunk of reads to process in parallel.
    #[arg(long, short='c', default_value_t = 300000)]
    parallel_chunk_size: usize,
    /// Edit distance used for reference sequences.
    #[arg(short, long, default_value_t = 20)]
    edit_dist: u8,
    /// FR4 region. Located just after the CDR3.
    #[arg(short, long, default_value_t = String::from("TGGGGCAAAGGGACCCAGGTCAC"))]
    fr4: String,
    /// Omit the header from the output. Useful if you're piping into another program to do more processing.
    #[arg(short='l', long, default_value_t = false)]
    headerless: bool
}

mod input {
    use super::*;

    pub enum Reader {
        BioFasta(Box<dyn BufRead>, usize),
        BioFastq(Box<dyn BufRead>, usize),
    }

    impl Reader {
        pub fn new(filename: &str, args: &Args) -> Result<Reader, InputError> {
            let (filetype, bufread) = reader(filename)?;

            Ok(match filetype {
                FileType::Fasta => Reader::BioFasta(bufread, args.parallel_chunk_size),
                FileType::Fastq => Reader::BioFastq(bufread, args.parallel_chunk_size),
            })
        }

        pub fn map<T: Send, R>(self, local_fn: impl Fn(Read) -> T + Sync, global_fn: impl Fn(T, &mut R), global_data: &mut R, head: Option<usize>) {
            match self {
                Reader::BioFasta(bufread, chunk_size) => {
                    let input_records = bio::io::fasta::Reader::from_bufread(bufread)
                        .records();

                    let recs = match head {
                        Some(n) => 
                            itertools::Either::Right(input_records.take(n)),
                        None => 
                            itertools::Either::Left(input_records)
                    };
        
                    for bunch in &recs.chunks(chunk_size) {
                        let mut outs = Vec::new();
                        
                        bunch.collect_vec().into_par_iter().map(|result| {
                            match result {
                                Err(_) => panic!("Bad record!"),
                                Ok(record) => local_fn(Read::BioFasta(&record))
                            }
                        }).collect_into_vec(&mut outs);
        
                        for t in outs {
                            global_fn(t, global_data)
                        }
                    }
                },

                Reader::BioFastq(bufread, chunk_size) => {
                    let input_records = bio::io::fastq::Reader::from_bufread(bufread)
                        .records();
        
                    for bunch in &input_records.chunks(chunk_size) {
                        let mut outs = Vec::new();
                        
                        bunch.collect_vec().into_par_iter().map(|result| {
                            match result {
                                Err(_) => panic!("Bad record!"),
                                Ok(record) => local_fn(Read::BioFastq(&record))
                            }
                        }).collect_into_vec(&mut outs);
        
                        for t in outs {
                            global_fn(t, global_data)
                        }
                    }
                },
            };
            

        }
    }

    use std::{ffi::OsStr, fs::File, io::{stdin, BufRead, BufReader}, path::Path, rc::Rc};
    use itertools::Itertools;

    use crate::Args;

    #[derive(Debug)]
    pub enum FileTypeError {
        UnknownFileType
    }

    pub enum FileType {
        Fasta,
        Fastq
    }

    impl TryFrom<&str> for FileType {
        type Error = FileTypeError;
    
        fn try_from(value: &str) -> Result<Self, Self::Error> {
            match value {
                "fa" => Ok(FileType::Fasta),
                "fasta" => Ok(FileType::Fasta),

                "fq" => Ok(FileType::Fastq),
                "fastq" => Ok(FileType::Fastq),

                _ => Err(FileTypeError::UnknownFileType),
            }
        }
    }

    #[derive(Debug)]
    pub enum InputError {
        FileNameError,
        FileTypeError(FileTypeError),
        FileOpenError(std::io::Error),
    }

    /// Checks the arguments, and either opens a file or reads from stdin
    pub fn reader(input_file: &str) -> Result<(FileType, Box<dyn BufRead>), InputError> {
        if input_file.eq("stdin") {
            // just read straight from stdin
            Ok((FileType::Fasta, Box::new(BufReader::new(stdin()))))
        } else {
            let path = Path::new(&input_file);

            let file = File::open(path)
                .map_err(InputError::FileOpenError)?;
            
            match &input_file.split('.').collect_vec()[..] {
                [.., ext, "gz"] =>
                    match FileType::try_from(*ext) {
                        Ok(filetype) => 
                            Ok((filetype, Box::new(
                                BufReader::new(flate2::read::MultiGzDecoder::new(file))))),
                        Err(err) => 
                            Err(InputError::FileTypeError(err)),
                    },
                
                [.., ext] =>
                    match FileType::try_from(*ext) {
                        Ok(filetype) => 
                            Ok((filetype, Box::new(
                                BufReader::new(file)))),
                        Err(err) => 
                            Err(InputError::FileTypeError(err)),
                    },
                
                _ => Err(InputError::FileNameError)
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

    pub fn print_header<T: Write>(output: &mut T) {
        writeln!(output, "id\tv_gene\tsequence\tcdr3_sequence\ttranslated_sequence")
            .expect("Couldn't write header line to output!");
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

enum Read<'a> {
    SeqIOFasta(seq_io::fasta::RefRecord<'a>),
    SeqIOFastq(seq_io::fastq::RefRecord<'a>),
    BioFasta(&'a bio::io::fasta::Record),
    BioFastq(&'a bio::io::fastq::Record),
}

use seq_io::fastq::Record as FastqRecord;
use seq_io::fasta::Record as FastaRecord;

impl <'a>Read<'a> {
    fn seq(&self) -> &[u8] {
        match self {
            Read::SeqIOFasta(fasta) => fasta.seq(),
            Read::SeqIOFastq(fastq) => fastq.seq(),

            Read::BioFasta(fasta) => fasta.seq(),
            Read::BioFastq(fastq) => fastq.seq(),
        }
    }

    fn name(&self) -> &str {
        match self {
            Read::SeqIOFasta(fasta) => fasta.id().ok().unwrap_or(""),
            Read::SeqIOFastq(fastq) => fastq.id().ok().unwrap_or(""),
            Read::BioFasta(fasta) => fasta.id(),
            Read::BioFastq(fastq) => fastq.id(),
        }
    } 
}

fn main() {
    bio_main();
}

// fn seq_io_main() -> ! {
//     // get the arguments from the command line
//     let args = Args::parse();

//     // collect and optimise the reference seqs
//     let a = input::reader(&args.input_reads);

//     let reference_seqs = reference::parse_reference(&args.reference_fasta);
//     let (optimised_reference_seqs, _rest) = reference::optimise_refs(
//         &reference_seqs, 
//         bio::io::fasta::Reader::from().records(),
//         args.sample_size, args.parallel_chunk_size, args.edit_dist, args.reference_size);
//     // let optimised_reference_seqs = reference_seqs;
    
//     let mut output_csv_file = output::writer(&args);
//     // print the header unless the user specifies not to
//     if !args.headerless {
//         output::print_header(&mut output_csv_file);
//     }

//     let fr4 = Myers::<u64>::new(args.fr4.as_bytes());

//     // read the input
//     let (filetype, reader) = input::reader(&args)
//         .unwrap();

//     match filetype {
//         input::FileType::Fasta => read_fasta(Box::new(reader), &optimised_reference_seqs[..], args, fr4, output_csv_file),
//         input::FileType::Fastq => (),
//     }
    
// }

// fn read_fastq() {

// }

// fn read_fasta(optimised_reference_seqs: &[RefV], args: Args, fr4: Myers, output_csv_file: Box<dyn Write>) {
    
//     let (a, b) = input::reader(&args).unwrap();
//     let reader = seq_io::fasta::Reader::(BufReader::new(flate2::read::MultiGzDecoder::new(File::create("hello.txt")))).unwrap();
    
//     let _ = parallel_fasta(reader, 64, 1000, 
//         |record, out| {
//             *out = find_cdr3::parse_read(
//                 &Read::Fasta(record),
//                 &optimised_reference_seqs, 
//                 args.edit_dist, 
//                 &fr4
//             );
//         }, 
//         |record, out| {
//             // either the seq, or the reverse complement of the seq
//             let seq = match out {
//                 find_cdr3::SmallOutputRecord::Reverse(_, _, _) => bio::alphabets::dna::revcomp(record.seq()),
//                 _ => Vec::from(record.seq()),
//             };
            
//             writeln!(output_csv_file, "{}\t{}\t{}",
//                 std::str::from_utf8(record.head()).unwrap(), 
//                 std::str::from_utf8(&seq).unwrap(), 
//                 out
//             ).expect("Couldn't write line to output!");

//             None::<()>
//         }
//     );
// }

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

fn bio_main() {
    // get the arguments from the command line
    let args = Args::parse();

    // collect all the reference seqs
    let reference_seqs = reference::parse_reference(&args.reference_fasta);

    // optimise the reference seqs list
    let (best_reference_seqs, _) = reference::optimise_refs(
        &reference_seqs,
        &args
    );

    let mut output_csv_file = output::writer(&args);
    // print the header unless the user specifies not to
    if !args.headerless {
        output::print_header(&mut output_csv_file);
    }

    let fr4 = Myers::<u64>::new(args.fr4.as_bytes());

    // go through each result and print to output while you go 
    // to avoid collecting data in memory
    let reader = input::Reader::new(&args.input_reads, &args)
        .unwrap();
    reader.map(
        |read| format!("{}\t{}\t{}", 
            String::from(read.name()),
            std::str::from_utf8(read.seq()).unwrap(),
            find_cdr3::parse_read(
                &read, 
                &best_reference_seqs,
                args.edit_dist, 
                &fr4
            )
        ),
        |out_string, output_csv_file| {
            writeln!(output_csv_file, "{out_string}")
                .expect("Couldn't write line to output!");
        },
        &mut output_csv_file,
        None
    );
}