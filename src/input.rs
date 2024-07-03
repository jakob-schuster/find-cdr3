use super::*;
use seq_io::fastq::Record as FastqRecord;
use seq_io::fasta::Record as FastaRecord;

/// A read from any accepted input file type.
/// Differently stored and accessed depending on parsing backend.
pub enum Read<'a> {
    SeqIOFasta(seq_io::fasta::RefRecord<'a>),
    SeqIOFastq(seq_io::fastq::RefRecord<'a>),

    BioFasta(&'a bio::io::fasta::Record),
    BioFastq(&'a bio::io::fastq::Record),
}


impl <'a>Read<'a> {
    /// Gets the seq of a read.
    pub fn seq(&self) -> &[u8] {
        match self {
            Read::SeqIOFasta(fasta) => fasta.seq(),
            Read::SeqIOFastq(fastq) => fastq.seq(),

            Read::BioFasta(fasta) => fasta.seq(),
            Read::BioFastq(fastq) => fastq.seq(),
        }
    }

    /// Gets the name of a read.
    pub fn name(&self) -> &str {
        match self {
            Read::SeqIOFasta(fasta) => fasta.id().ok().unwrap_or(""),
            Read::SeqIOFastq(fastq) => fastq.id().ok().unwrap_or(""),
            Read::BioFasta(fasta) => fasta.id(),
            Read::BioFastq(fastq) => fastq.id(),
        }
    } 
}

/// A reader for input files, taking either Fastq or Fasta files.
/// Currently, only the rust-bio parser backend is used.
pub enum Reader {
    BioFasta(Box<dyn BufRead>, usize, Option<usize>),
    BioFastq(Box<dyn BufRead>, usize, Option<usize>),
}

impl Reader {
    /// Creates a new reader, given the filename, 
    /// the size of each parallel-processed chunk of reads, 
    /// and the option to only take the first n reads.
    pub fn new(
        filename: &str, 
        parallel_chunk_size: usize, 
        take_first: Option<usize>
    ) -> Result<Reader, InputError> {
        let (filetype, bufread) = reader(filename)?;

        // just use bio parsers for now; can't use seq_io in parallel on gzipped files
        Ok(match filetype {
            FileType::Fasta => Reader::BioFasta(bufread, parallel_chunk_size, take_first),
            FileType::Fastq => Reader::BioFastq(bufread, parallel_chunk_size, take_first),
        })
    }

    /// Maps a function across all the reads in the input.
    pub fn map<T: Send, R>(self, local_fn: impl Fn(Read) -> T + Sync, global_fn: impl Fn(T, &mut R), global_data: &mut R) {
        match self {
            Reader::BioFasta(bufread, chunk_size, take_first) => {
                let input_records = bio::io::fasta::Reader::from_bufread(bufread)
                    .records();

                let recs = match take_first {
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

            Reader::BioFastq(bufread, chunk_size, take_first) => {
                let input_records = bio::io::fastq::Reader::from_bufread(bufread)
                    .records();
    
                let recs = match take_first {
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

/// Checks the arguments, and either opens a file or reads from stdin.
/// Also returns the filetype.
fn reader(input_file: &str) -> Result<(FileType, Box<dyn BufRead>), InputError> {
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