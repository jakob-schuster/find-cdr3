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

fn find_cdr3(seq: &String, reference_seqs: &Vec<reference::RefV>, edit_dist: u8) -> Option<Vec<u8>> {
    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_seqs
    .into_iter().map(|reference::RefV { seq: v_seq, myers, cys_index } | {
        let mut m2 = myers.to_owned();

        let matches = m2.find_all_lazy(seq.as_bytes(), 8);

        let best_match = matches.into_iter().map(|(end_pos, _)| {
            let mut aln = Alignment::default();
            // let mut ops: Vec<_> = Vec::new();
            matches.alignment_at(end_pos, &mut aln);
            (aln, edit_dist)
        }).max_by_key(|(_, edit_dist)| *edit_dist);

        let a = match best_match {
            None => None,
            Some((aln, b)) => {
                let (cys_start_seq, cys_start_pat) = aln.to_owned().operations.into_iter()
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
                                AlignmentOperation::Xclip(n) => {panic!("hit an xclip"); (seq_i + n, v_i)},
                                AlignmentOperation::Yclip(n) => {panic!("hit a yclip"); (seq_i, v_i + n)}
                            }
                        }
                    }
                );
                assert_eq!(cys_start_pat, *cys_index);

                // return the calculated cys start
                Some((cys_start_seq + aln.ystart, edit_dist))
            } 

            // .map_or(None, |(a, b)| {
            //     let mut ops: Vec<_> = Vec::new(); 
            //     matches.path_at(a, &mut ops); 
                
            //     let (cys_start_seq, cys_start_pat) = aln.operations.into_iter()
            //         .fold((0, 0), |(seq_i, v_i), operation| {
            //             if v_i == *cys_index {
            //                 // already reached the cys; just pass it through unchanged
            //                 (seq_i, v_i)
            //             } else {
            //                 // handle the current operation
            //                 match operation {
            //                     AlignmentOperation::Del => (seq_i, v_i + 1),
            //                     AlignmentOperation::Ins => (seq_i + 1, v_i),
            //                     AlignmentOperation::Match => (seq_i + 1, v_i + 1),
            //                     AlignmentOperation::Subst => (seq_i + 1, v_i + 1),
            //                     AlignmentOperation::Xclip(n) => {panic!("hit an xclip"); (seq_i + n, v_i)},
            //                     AlignmentOperation::Yclip(n) => {panic!("hit a yclip"); (seq_i, v_i + n)}
            //                 }
            //             }
            //         }
            //     );
            //     assert_eq!(cys_start_pat, *cys_index);

            //     // return the calculated cys start
            //     Some((cys_start_seq + aln.ystart, edit_dist))
            // });

            best_match
        }

        // let best_match = myers.find_best_end(seq.as_bytes());


        // let best_match = matches.by_ref().max_by_key(|&(_, _, dist)| dist);

        // match best_match {
        //     None => None,
        //     Some((_, best_start, best_end, best_dist)) => {
        //         matches.path_at(best_start);
        //         // the start of the cys codon in the actual sequence, 
        //         // according to the position of cys in the reference 
        //         // and the operations involved in the alignment

        //         let ops = aln.operations.clone();
        //         let cys_start_seq = aln.operations.into_iter()
        //             .fold((0, 0), |(seq_i, v_i), operation| {
        //                 if v_i == *cys_index {
        //                     // already reached the cys; just pass it through unchanged
        //                     (seq_i, v_i)
        //                 } else {
        //                     // handle the current operation
        //                     match operation {
        //                         AlignmentOperation::Del => (seq_i, v_i + 1),
        //                         AlignmentOperation::Ins => (seq_i + 1, v_i),
        //                         AlignmentOperation::Match => (seq_i + 1, v_i + 1),
        //                         AlignmentOperation::Subst => (seq_i + 1, v_i + 1),
        //                         AlignmentOperation::Xclip(n) => {panic!("hit an xclip"); (seq_i + n, v_i)},
        //                         AlignmentOperation::Yclip(n) => {panic!("hit a yclip"); (seq_i, v_i + n)}
        //                     }
        //                 }
        //             }
        //         ).0;
                
        //         let start = aln.ystart;
        // Some((cys_start_seq + aln.ystart, edit_dist))
            
        //     }
        // }
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
            match fr4_imgt_myers.find_all(&seq.as_bytes()[cys_start_seq..], edit_dist / 2)
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