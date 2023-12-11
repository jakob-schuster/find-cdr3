use std;

use std::fs::File;

use std::io::BufReader;


use bio;

use std::cmp;

use bio::alignment::AlignmentOperation;

use bio::alignment::Alignment;

use bio::pattern_matching::myers::Myers;

use core::fmt;

use crate::reference;

pub(crate) struct OutputRecord {
    pub(crate) name: String,
    pub(crate) seq: String,
    pub(crate) cdr3_seq: String
}

impl fmt::Display for OutputRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.name, self.seq, self.cdr3_seq)
    }
}

pub(crate) fn find_cdr3(seq: &String, reference_seqs: &Vec<reference::RefV>, edit_dist: u8) -> Option<Vec<u8>> {
    // build ALL the myers matchers inside each iteration - bad?
    let reference_myers: Vec<(Myers, &usize)> = reference_seqs.into_iter()
        .map(|reference::RefV { seq, cys_index } | 
            (Myers::<u64>::new(seq.as_str().as_bytes()), cys_index))
        .collect();

    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_myers
    .into_iter().map(|(mut myers, cys_index)| {
        let mut aln = Alignment::default();
    
        let mut matches = myers.find_all_lazy(seq.as_bytes(), edit_dist);
        let best_match = matches.by_ref().max_by_key(|&(_, dist)| dist);

        match best_match {
            None => None,
            Some((best_end, best_dist)) => {
                matches.alignment_at(best_end, &mut aln);
            
                // the start of the cys codon in the actual sequence, 
                // according to the position of cys in the reference 
                // and the operations involved in the alignment
                let cys_start_seq = aln.operations.into_iter()
                    .fold((0, 0), |(seq_i, v_i), operation| {
                        if v_i == *cys_index {
                            // already reached the cys; just pass it through unchanged
                            (seq_i, v_i)
                        } else {
                            // handle the current operation
                            match operation {
                                AlignmentOperation::Del => (seq_i, v_i + 1),
                                AlignmentOperation::Ins => (seq_i + 1, v_i),
                                AlignmentOperation::Match => (seq_i + 1, v_i + 1),
                                AlignmentOperation::Subst => (seq_i + 1, v_i + 1),
                                AlignmentOperation::Xclip(n) => (seq_i + n, v_i),
                                AlignmentOperation::Yclip(n) => (seq_i, v_i + n)
                            }
                        }
                    }
                ).0;

                let start = aln.ystart;
            
                Some((cys_start_seq + aln.ystart, best_dist))
            }
        }
    });


    // get the cys start of the best match, if there is a match at all
    let v_best_match = v_matches.max_by(|a, b| {
        match (a, b) {
            (None, None) => cmp::Ordering::Equal,
            (None, Some(_)) => cmp::Ordering::Less,
            (Some(_), None) => cmp::Ordering::Greater,
            (Some((_, edit_dist_a)), Some((_, edit_dist_b))) => edit_dist_a.cmp(edit_dist_b) 
        }
    }).map_or(None, |a| {
        a.map_or(None, |(cys_start, _)| {
            Some(cys_start)
        })
    });

    // currently hardcoded codon for the fr4
    let mut fr4_imgt_myers = Myers::<u64>::new(b"TGGGGCAAAGGGACCCAGGTCAC");

    match v_best_match {
        None => None,
        Some(cys_start_seq) => {
            match fr4_imgt_myers.find_all(&seq.as_bytes()[cys_start_seq..], 3)
            .min_by_key(|&(start, _, _)| start) {
                None => None,
                Some(first_fr4_match) => {
                    let fr4_start = first_fr4_match.0 + cys_start_seq;
                    Some(seq.as_bytes()[(cys_start_seq+3)..fr4_start].to_owned())
                }
            }
        } 
    }
}

pub(crate) fn parse_input(input_fasta: String, reference_seqs: &Vec<reference::RefV>, edit_dist: u8) -> Vec<OutputRecord> {
    let reader = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(input_fasta).expect("Couldn't open input!")));

    let recs: Vec<_> = reader.records().collect();

    recs.into_iter()
        .map(|result| match result {
            Err(_) => panic!("Bad record in input!"),
            Ok(record) => {
                let forward = find_cdr3(&String::from_utf8(
                    record.seq().to_owned()).unwrap(), 
                    &reference_seqs, 
                    edit_dist
                );
                let reverse = find_cdr3(
                    &String::from_utf8(bio::alphabets::dna::revcomp(
                        record.seq().to_owned())).unwrap(),
                    &reference_seqs, 
                    edit_dist
                );

                if forward > reverse {
                    OutputRecord { 
                        name: record.id().to_string(),
                        seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                        cdr3_seq: String::from_utf8(forward.unwrap()).unwrap()
                    }
                } else if reverse > forward {
                    OutputRecord { 
                        name: record.id().to_string(),
                        seq: String::from_utf8(bio::alphabets::dna::revcomp(
                            record.seq().to_owned())).unwrap(),
                        cdr3_seq: String::from_utf8(reverse.unwrap()).unwrap()
                    }
                } else {
                    OutputRecord {
                        name: record.id().to_string(),
                        seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                        cdr3_seq: "".to_string()
                    }
                }
            }
        }
    ).collect()
}