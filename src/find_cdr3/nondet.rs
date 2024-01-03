use std::{cmp, thread};

use bio::{alignment::{Alignment, AlignmentOperation}, pattern_matching::myers::Myers, io::fasta::Record};
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IndexedParallelIterator};

use crate::reference;

use super::OutputRecord;

pub(crate) fn parse_one_input(
    record: Record,
    reference_seqs: &Vec<reference::RefV>,
    edit_dist: u8
) -> OutputRecord {
    let forward = find_cdr3(
        record.seq(), 
        &reference_seqs,
        edit_dist
    );

    let reverse = find_cdr3(
        &bio::alphabets::dna::revcomp(record.seq()),
        &reference_seqs,
        edit_dist
    );

    match (&forward, &reverse) {
        // neither was successful
        (None, None) => OutputRecord {
            name: String::from(record.id()),
            seq: String::from_utf8(record.seq().to_vec()).unwrap(),
            cdr3_seq: String::new()
        },

        // forwards was successful
        (Some(seq), None) => OutputRecord { 
            name: String::from(record.id()),
            seq: String::from_utf8(record.seq().to_vec()).unwrap(),
            cdr3_seq: String::from_utf8(seq.to_owned()).unwrap()
        },

        // reverse was successful
        (None, Some(seq)) => OutputRecord { 
            name: String::from(record.id()),
            seq: String::from_utf8(bio::alphabets::dna::revcomp(record.seq())).unwrap(),
            cdr3_seq: String::from_utf8(seq.to_owned()).unwrap()
        },

        // both were successful - ambiguous
        (Some(seq), Some(rev_seq)) => OutputRecord {
            name: String::from(record.id()),
            seq: String::from_utf8(record.seq().to_vec()).unwrap(),
            cdr3_seq: String::new()
        }
    }
}

pub(crate) fn parse_one_input_par(
    record: Record,
    reference_seqs: &Vec<reference::RefV>,
    edit_dist: u8
) -> OutputRecord {
    let mut vec = Vec::new();
    
    vec![record.seq(), &bio::alphabets::dna::revcomp(record.seq())]
        .into_par_iter()
        .map(|seq| {
            find_cdr3_par(
                seq, 
                &reference_seqs.to_owned(),
                edit_dist
            )
        })
        .collect_into_vec(&mut vec);

    if let [forward, reverse] = &vec[..] {
        match (&forward, &reverse) {
            // neither was successful
            (None, None) => OutputRecord {
                name: String::from(record.id()),
                seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                cdr3_seq: String::new()
            },
    
            // forwards was successful
            (Some(seq), None) => OutputRecord { 
                name: String::from(record.id()),
                seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                cdr3_seq: String::from_utf8(seq.to_owned()).unwrap()
            },
    
            // reverse was successful
            (None, Some(seq)) => OutputRecord { 
                name: String::from(record.id()),
                seq: String::from_utf8(bio::alphabets::dna::revcomp(record.seq())).unwrap(),
                cdr3_seq: String::from_utf8(seq.to_owned()).unwrap()
            },
    
            // both were successful - ambiguous
            (Some(seq), Some(rev_seq)) => OutputRecord {
                name: String::from(record.id()),
                seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                cdr3_seq: String::new()
            }
        }
    } else {
        panic!("Thread problem")
    }
}

pub(crate) fn find_cdr3(
    seq: &[u8], reference_seqs: &Vec<reference::RefV>, edit_dist: u8
) -> Option<Vec<u8>> {
    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_seqs.into_iter()
    .map(|reference::RefV { seq: _, myers, cys_index } | {
        let mut aln = Alignment::default();
    
        let mut owned_myers = myers.clone();

        let mut matches = 
            owned_myers.find_all_lazy(seq, edit_dist);
        let best_match = 
            matches.by_ref().max_by_key(|&(_, dist)| dist);

        match best_match {
            None => None,
            Some((best_end, best_dist)) => {
                matches.alignment_at(best_end, &mut aln);

                // the end of the cys index within the reference
                let cys_index_end = *cys_index + 3;
                // the start of the cys codon in the actual sequence, 
                // according to the position of cys in the reference 
                // and the operations involved in the alignment
                let (mut cys_end_seq, mut cys_start_seq):
                    (Option<usize>, Option<usize>) = (None, None);
                let (mut seq_i, mut v_i): (usize, usize) = (0, 0);

                for op in aln.operations {
                    if v_i == *cys_index {
                        cys_start_seq = Some(seq_i);
                    } else if v_i == cys_index_end {
                        cys_end_seq = Some(seq_i);
                        break;
                    }

                    match op {
                        AlignmentOperation::Del => { v_i += 1 },
                        AlignmentOperation::Ins => { seq_i += 1 },
                        AlignmentOperation::Match => { seq_i += 1; v_i += 1 },
                        AlignmentOperation::Subst => { seq_i += 1; v_i += 1 },
                        AlignmentOperation::Xclip(_) => panic!("Xclip"),
                        AlignmentOperation::Yclip(_) => panic!("Yclip")
                    }
                }
                
                match (cys_start_seq, cys_end_seq) {
                    (Some(n), Some(_)) =>
                        // Some((n + aln.ystart, m + aln.ystart, best_dist)),
                        Some((n + aln.ystart, n + 3 + aln.ystart, best_dist)), // cheap hack
                    _ => 
                        None
                }
            }
        }
    });

    // get the cys start of the best match, if there is a match at all
    let v_best_match = v_matches.max_by(|a, b| {
        match (a, b) {
            (None, None) => cmp::Ordering::Equal,
            (None, Some(_)) => cmp::Ordering::Less,
            (Some(_), None) => cmp::Ordering::Greater,
            (Some((_, _, edit_dist_a)), Some((_, _, edit_dist_b))) 
                => edit_dist_a.cmp(edit_dist_b)
        }
    }).map_or(None, |a| {
        a.map_or(None, |(cys_start, cys_end, _)| {
            Some((cys_start, cys_end))
        })
    });

    // currently hardcoded codon for the fr4
    let mut fr4_imgt_myers = Myers::<u64>::new(b"TGGGGCAAAGGGACCCAGGTCAC");

    match v_best_match {
        None => None,
        Some((cys_start_seq, cys_end_seq)) => {
            if cys_end_seq >= seq.len() {
                // some reads get cut off just after the V gene (that is, FR3-IMGT)
                // println!("Seq {} len {}. Cys at ({},{}).", 
                //     seq.to_ascii_uppercase(), seq.len(), cys_start_seq, cys_end_seq);
                None
            } else {
                match fr4_imgt_myers
                    .find_all(&seq[cys_end_seq..], 3)
                    .min_by_key(|&(start, _, _)| start) {
                        None => None,
                        Some(first_fr4_match) => {
                            // println!("Seq {} len {}. Cys at ({},{}).", seq, seq.len(), cys_start_seq, cys_end_seq);
                            let fr4_start = first_fr4_match.0 + cys_end_seq;
                            Some(seq[cys_end_seq..fr4_start]
                                .to_owned())
                        }
                    }
            }
        } 
    }
}

pub(crate) fn find_cdr3_par(
    seq: &[u8], reference_seqs: &Vec<reference::RefV>, edit_dist: u8
) -> Option<Vec<u8>> {
    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_seqs.into_par_iter()
    .map(|reference::RefV { seq: _, myers, cys_index } | {
        let mut aln = Alignment::default();
    
        let mut owned_myers = myers.clone();

        let mut matches = 
            owned_myers.find_all_lazy(seq, edit_dist);
        let best_match = 
            matches.by_ref().max_by_key(|&(_, dist)| dist);

        match best_match {
            None => None,
            Some((best_end, best_dist)) => {
                matches.alignment_at(best_end, &mut aln);

                // the end of the cys index within the reference
                let cys_index_end = *cys_index + 3;
                // the start of the cys codon in the actual sequence, 
                // according to the position of cys in the reference 
                // and the operations involved in the alignment
                let (mut cys_end_seq, mut cys_start_seq):
                    (Option<usize>, Option<usize>) = (None, None);
                let (mut seq_i, mut v_i): (usize, usize) = (0, 0);

                for op in aln.operations {
                    if v_i == *cys_index {
                        cys_start_seq = Some(seq_i);
                    } else if v_i == cys_index_end {
                        cys_end_seq = Some(seq_i);
                        break;
                    }

                    match op {
                        AlignmentOperation::Del => { v_i += 1 },
                        AlignmentOperation::Ins => { seq_i += 1 },
                        AlignmentOperation::Match => { seq_i += 1; v_i += 1 },
                        AlignmentOperation::Subst => { seq_i += 1; v_i += 1 },
                        AlignmentOperation::Xclip(_) => panic!("Xclip"),
                        AlignmentOperation::Yclip(_) => panic!("Yclip")
                    }
                }
                
                match (cys_start_seq, cys_end_seq) {
                    (Some(n), Some(_)) =>
                        // Some((n + aln.ystart, m + aln.ystart, best_dist)),
                        Some((n + aln.ystart, n + 3 + aln.ystart, best_dist)), // cheap hack
                    _ => 
                        None
                }
            }
        }
    });

    // get the cys start of the best match, if there is a match at all
    let v_best_match = v_matches.max_by(|a, b| {
        match (a, b) {
            (None, None) => cmp::Ordering::Equal,
            (None, Some(_)) => cmp::Ordering::Less,
            (Some(_), None) => cmp::Ordering::Greater,
            (Some((_, _, edit_dist_a)), Some((_, _, edit_dist_b))) 
                => edit_dist_a.cmp(edit_dist_b)
        }
    }).map_or(None, |a| {
        a.map_or(None, |(cys_start, cys_end, _)| {
            Some((cys_start, cys_end))
        })
    });

    // currently hardcoded codon for the fr4
    let mut fr4_imgt_myers = Myers::<u64>::new(b"TGGGGCAAAGGGACCCAGGTCAC");

    match v_best_match {
        None => None,
        Some((cys_start_seq, cys_end_seq)) => {
            if cys_end_seq >= seq.len() {
                // some reads get cut off just after the V gene (that is, FR3-IMGT)
                // println!("Seq {} len {}. Cys at ({},{}).", 
                //     seq.to_ascii_uppercase(), seq.len(), cys_start_seq, cys_end_seq);
                None
            } else {
                match fr4_imgt_myers
                    .find_all(&seq[cys_end_seq..], 3)
                    .min_by_key(|&(start, _, _)| start) {
                        None => None,
                        Some(first_fr4_match) => {
                            // println!("Seq {} len {}. Cys at ({},{}).", seq, seq.len(), cys_start_seq, cys_end_seq);
                            let fr4_start = first_fr4_match.0 + cys_end_seq;
                            Some(seq[cys_end_seq..fr4_start]
                                .to_owned())
                        }
                    }
            }
        } 
    }
}
