use super::*;
use seq_io::fasta::Record as FastaRecord;
use seq_io::fastq::Record as FastqRecord;

/// A read from any accepted input file type.
/// Differently stored and accessed depending on parsing backend.
pub enum Read<'a> {
    SeqIOFasta(seq_io::fasta::RefRecord<'a>),
    SeqIOFastq(seq_io::fastq::RefRecord<'a>),

    BioFasta(&'a bio::io::fasta::Record),
    BioFastq(&'a bio::io::fastq::Record),
}

impl<'a> Read<'a> {
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
    BioFasta {
        bufread: Box<dyn BufRead>,
        chunk_size: usize,
        threads: usize,
        take_first: Option<usize>,
    },
    BioFastq {
        bufread: Box<dyn BufRead>,
        chunk_size: usize,
        threads: usize,
        take_first: Option<usize>,
    },
    NonGzippedSeqIOFastq {
        filename: String,
        threads: usize,
    },
    NonGzippedSeqIOFasta {
        filename: String,
        threads: usize,
    },
}

impl Reader {
    pub fn new(
        filename: &str,
        chunk_size: usize,
        threads: usize,
        take_first: Option<usize>,
    ) -> Result<Reader, InputError> {
        Self::new_bio(filename, chunk_size, threads, take_first)
    }

    /// Creates a new reader, given the filename,
    /// the size of each parallel-processed chunk of reads,
    /// and the option to only take the first n reads.
    /// Always just uses the bio reader.
    fn new_bio(
        filename: &str,
        chunk_size: usize,
        threads: usize,
        take_first: Option<usize>,
    ) -> Result<Reader, InputError> {
        let (filetype, bufread) = reader(filename)?;

        Ok(match filetype.inner() {
            RawFileType::Fasta => Reader::BioFasta {
                bufread,
                chunk_size,
                threads,
                take_first,
            },
            RawFileType::Fastq => Reader::BioFastq {
                bufread,
                chunk_size,
                threads,
                take_first,
            },
        })
    }

    /// Creates a new reader, given the filename,
    /// the size of each parallel-processed chunk of reads,
    /// and the option to only take the first n reads.
    /// Chooses between a seq_io reader, and a bio reader,
    /// based on whether the file is gzipped or not
    /// (gzipped parallel processing is not compatible with
    /// seq_io). This was ditched because it turns out there's no real performance difference.
    fn new_choose(
        filename: &str,
        chunk_size: usize,
        threads: usize,
        take_first: Option<usize>,
    ) -> Result<Reader, InputError> {
        let (filetype, bufread) = reader(filename)?;

        Ok(match take_first {
            // if we're taking just the first reads, always use bio
            Some(_) => Self::new_bio(filename, chunk_size, threads, take_first)?,
            None => {
                // otherwise, choose intelligently
                match filetype {
                    // if not gzipped, we can use seq_io
                    FileType::Raw(raw) => match raw {
                        RawFileType::Fasta => Reader::NonGzippedSeqIOFasta {
                            filename: filename.to_string(),
                            threads,
                        },
                        RawFileType::Fastq => Reader::NonGzippedSeqIOFastq {
                            filename: filename.to_string(),
                            threads,
                        },
                    },

                    // otherwise, have to use bio
                    FileType::Gzipped(_) => {
                        Self::new_bio(filename, chunk_size, threads, take_first)?
                    }
                }
            }
        })
    }

    /// Maps a function across all the reads in the input.
    pub fn map<T: Send + Default, R>(
        self,
        local_fn: impl Fn(Read) -> T + Sync,
        global_fn: impl Fn(&T, &mut R),
        global_data: &mut R,
    ) {
        match self {
            Reader::BioFasta {
                bufread,
                chunk_size,
                threads,
                take_first,
            } => {
                let input_records = bio::io::fasta::Reader::from_bufread(bufread).records();

                let recs = match take_first {
                    Some(n) => itertools::Either::Right(input_records.take(n)),
                    None => itertools::Either::Left(input_records),
                };

                for bunch in &recs.chunks(chunk_size) {
                    let mut outs = Vec::new();

                    bunch
                        .collect_vec()
                        .into_par_iter()
                        .map(|result| match result {
                            Err(_) => panic!("Bad record!"),
                            Ok(record) => local_fn(Read::BioFasta(&record)),
                        })
                        .collect_into_vec(&mut outs);

                    for t in outs {
                        global_fn(&t, global_data)
                    }
                }
            }

            Reader::BioFastq {
                bufread,
                chunk_size,
                threads,
                take_first,
            } => {
                let input_records = bio::io::fastq::Reader::from_bufread(bufread).records();

                let recs = match take_first {
                    Some(n) => itertools::Either::Right(input_records.take(n)),
                    None => itertools::Either::Left(input_records),
                };

                for bunch in &recs.chunks(chunk_size) {
                    let mut outs = Vec::new();

                    bunch
                        .collect_vec()
                        .into_par_iter()
                        .map(|result| match result {
                            Err(_) => panic!("Bad record!"),
                            Ok(record) => local_fn(Read::BioFastq(&record)),
                        })
                        .collect_into_vec(&mut outs);

                    for t in outs {
                        global_fn(&t, global_data)
                    }
                }
            }

            Reader::NonGzippedSeqIOFasta { filename, threads } => {
                let input_records = seq_io::fasta::Reader::from_path(filename).unwrap();

                seq_io::parallel::parallel_fasta(
                    input_records,
                    threads as u32,
                    10,
                    |rec, out| *out = local_fn(Read::SeqIOFasta(rec)),
                    |_, out| {
                        global_fn(out, global_data);
                        Some(())
                    },
                )
                .unwrap();
            }

            Reader::NonGzippedSeqIOFastq { filename, threads } => {
                let input_records = seq_io::fastq::Reader::from_path(filename).unwrap();

                seq_io::parallel::parallel_fastq(
                    input_records,
                    threads as u32,
                    10,
                    |rec, out| *out = local_fn(Read::SeqIOFastq(rec)),
                    |_, out| {
                        global_fn(out, global_data);
                        Some(())
                    },
                )
                .unwrap();
            }
        };
    }
}

use itertools::Itertools;
use std::{
    fs::File,
    io::{stdin, BufRead, BufReader},
    path::Path,
};

#[derive(Debug)]
pub enum FileTypeError {
    UnknownFileType,
}

pub enum FileType {
    Raw(RawFileType),
    Gzipped(RawFileType),
}
impl FileType {
    fn inner(&self) -> &RawFileType {
        match self {
            FileType::Raw(r) => r,
            FileType::Gzipped(r) => r,
        }
    }
}

pub enum RawFileType {
    Fasta,
    Fastq,
}

impl TryFrom<&str> for RawFileType {
    type Error = FileTypeError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "fa" => Ok(RawFileType::Fasta),
            "fasta" => Ok(RawFileType::Fasta),

            "fq" => Ok(RawFileType::Fastq),
            "fastq" => Ok(RawFileType::Fastq),

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
        Ok((
            FileType::Raw(RawFileType::Fasta),
            Box::new(BufReader::new(stdin())),
        ))
    } else {
        let path = Path::new(&input_file);

        let file = File::open(path).map_err(InputError::FileOpenError)?;

        match &input_file.split('.').collect_vec()[..] {
            [.., ext, "gz"] => match RawFileType::try_from(*ext) {
                Ok(filetype) => Ok((
                    FileType::Gzipped(filetype),
                    Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))),
                )),
                Err(err) => Err(InputError::FileTypeError(err)),
            },

            [.., ext] => match RawFileType::try_from(*ext) {
                Ok(filetype) => Ok((FileType::Gzipped(filetype), Box::new(BufReader::new(file)))),
                Err(err) => Err(InputError::FileTypeError(err)),
            },

            _ => Err(InputError::FileNameError),
        }
    }
}
