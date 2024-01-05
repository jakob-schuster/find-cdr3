use std;

use std::fs::File;

use std::io::BufReader;
use std::io::Read;

use bio;
use bio::io::fasta::Record;

use std::cmp;

use bio::alignment::AlignmentOperation;

use bio::alignment::Alignment;

use bio::pattern_matching::myers::Myers;

use core::fmt;

use crate::reference;

pub struct OutputRecord {
    pub(crate) name: String,
    pub(crate) seq: String,
    pub(crate) cdr3_seq: String
}

impl fmt::Display for OutputRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.name, self.seq, self.cdr3_seq)
    }
}

pub mod det;
pub mod nondet;

pub(crate) fn parse_input(
    input_fasta: &str, 
    reference_seqs: &Vec<reference::RefV>, 
    edit_dist: u8,
    fr4: &[u8]
) -> Vec<OutputRecord> {
    let reader = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(input_fasta).expect("Couldn't open input!")));

    reader.records().map(|result| match result {
            Err(_) => panic!("Bad record in input!"),
            Ok(record) => {
                let forward = det::find_cdr3(
                    record.seq(), 
                    &reference_seqs, 
                    edit_dist,
                    fr4
                );
                
                if let Some(seq) = forward {
                    OutputRecord {
                        name: String::from(record.id()),
                        seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                        cdr3_seq: String::from_utf8(seq).unwrap()
                    }
                } else {
                    let reverse = det::find_cdr3(
                        &bio::alphabets::dna::revcomp(record.seq()),
                        &reference_seqs, 
                        edit_dist,
                        fr4
                    );

                    if let Some(seq) = reverse {
                        OutputRecord { 
                            name: String::from(record.id()),
                            seq: String::from_utf8(bio::alphabets::dna::revcomp(record.seq())).unwrap(),
                            cdr3_seq: String::from_utf8(seq).unwrap()
                        }
                    } else {
                        OutputRecord {
                            name: String::from(record.id()),
                            seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                            cdr3_seq: String::new()
                        }
                    }
                }
            }
        }
    ).collect()
}