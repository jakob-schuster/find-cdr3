use core::fmt;
use std::{cmp, fmt::Display};

use bio::{pattern_matching::myers::Myers, alignment::{AlignmentOperation, Alignment}};
use itertools::Itertools;

use crate::reference;

#[derive(Default)]
pub struct OutputRecord {
    pub(crate) name: String,
    pub(crate) seq: String,
    pub(crate) cdr3_seq: String
}

impl fmt::Display for OutputRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.name, self.seq, self.cdr3_seq)
    }
}

#[derive(Default, Debug)]
pub enum SmallOutputRecord {
    #[default]
    Neither,
    Forward(String, Vec<u8>, Vec<u8>),
    Reverse(String, Vec<u8>, Vec<u8>),
    Both
}

impl Display for SmallOutputRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SmallOutputRecord::Neither => 
                write!(f, "\t\t\t"),
            SmallOutputRecord::Forward(name, cdr3, full_seq) =>
                write!(f, "{}\t{}\t{}",
                    name,
                    std::str::from_utf8(cdr3).unwrap(),
                    // std::str::from_utf8(full_seq).unwrap(),
                    crate::reference::translate(full_seq)
                ),
            SmallOutputRecord::Reverse(name, cdr3, full_seq) =>
                write!(f, "{}\t{}\t{}",
                    name,
                    std::str::from_utf8(cdr3).unwrap(),
                    // std::str::from_utf8(full_seq).unwrap(),
                    crate::reference::translate(full_seq)
                ),
            SmallOutputRecord::Both =>
                write!(f, "\t\t\t"),
        }
    }
}

use super::*;

pub fn parse_read(
    read: &Read,
    reference_seqs: &[reference::RefV],
    edit_dist: u8,
    fr4: &Myers::<u64>
) -> SmallOutputRecord {
    let seq = read.seq();
    let rev_seq = &bio::alphabets::dna::revcomp(seq);

    let vec = [seq, rev_seq]
        .iter()
        .map(|seq| {
            find_cdr3(seq, reference_seqs, edit_dist, fr4)
        }).collect_vec();

    match &vec[..] {
        // neither was successful
        [None, None] => SmallOutputRecord::Neither,

        // forwards was successful
        [Some((name, seq_cdr3, full_seq)), None] => 
            SmallOutputRecord::Forward(name.clone(), seq_cdr3.clone(), full_seq.clone()),

        // reverse was successful
        [None, Some((name, rev_seq_cdr3, rev_full_seq))] => 
            SmallOutputRecord::Reverse(name.clone(), rev_seq_cdr3.clone(), rev_full_seq.clone()),

        // both were successful - ambiguous
        [Some(_), Some(_)] => SmallOutputRecord::Both,

        _ => panic!("Thread problem")
    }
}

fn find_cdr3(
    seq: &[u8], reference_seqs: &[reference::RefV], edit_dist: u8, fr4: &Myers::<u64>
) -> Option<(String, Vec<u8>, Vec<u8>)> {
    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_seqs.iter()
        .map(|reference::RefV { name, seq: _ref_seq, myers, cys_index: cys_end_ref } | {
            let mut owned_myers = myers.clone();

            let mut matches = 
                owned_myers.find_all_lazy(seq, edit_dist);
            let (best_end, _) = 
                matches.by_ref().min_by_key(|&(_, dist)| dist)?;

            let mut aln = Alignment::default();
            matches.alignment_at(best_end, &mut aln);

            Some((name, aln, cys_end_ref))
        });

    // get the cys start of the best match, if there is a match at all
    let (v_name, best_aln, cys_end_ref) = v_matches.min_by(|a, b| {
        match (a, b) {
            (None, None) => cmp::Ordering::Equal,
            (None, Some(_)) => cmp::Ordering::Greater,
            (Some(_), None) => cmp::Ordering::Less,
            (Some((_, aln_a, _)), Some((_, aln_b, _))) =>
                aln_a.score.cmp(&aln_b.score)
        }
    })??;

    // previously, we've calculated this properly, using the alignment path
    // igblast instead just finds the start position, then goes 3 bases further
    // so we've mimicked their technique
    let cys_start_seq = find_index(&best_aln, &(cys_end_ref-3))?;
    let cys_end_seq = cys_start_seq + 3;

    if cys_end_seq >= seq.len() {
        // some reads get cut off just after the V gene (that is, FR3-IMGT)
        None
    } else {
        let (best_fr4_start, best_fr4_end, _) = fr4.clone()
            .find_all(&seq[cys_end_seq..], edit_dist)
            .min_by_key(|&(_, _, dist)| dist)?;
        let fr4_start = best_fr4_start + cys_end_seq;
        
        Some((v_name.clone(), Vec::from(&seq[cys_end_seq..fr4_start]), Vec::from(&seq[best_aln.ystart..best_fr4_end + cys_end_seq])))
    }
}


fn find_cdr3_no_clone<'a>(
    seq: &'a [u8], 
    reference_seqs: &[reference::RefV], 
    edit_dist: u8, 
    fr4: &Myers::<u64>
) -> Option<&'a[u8]> {
    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_seqs.iter()
        .map(|reference::RefV { name, seq: _ref_seq, myers, cys_index: cys_end_ref } | {
            let mut owned_myers = myers.clone();

            let mut matches = 
                owned_myers.find_all_lazy(seq, edit_dist);
            let (best_end, _) = 
                matches.by_ref().min_by_key(|&(_, dist)| dist)?;

            let mut aln = Alignment::default();
            matches.alignment_at(best_end, &mut aln);

            Some((name, aln, cys_end_ref))
        });

    // get the cys start of the best match, if there is a match at all
    let (_, best_aln, cys_end_ref) = v_matches.min_by(|a, b| {
        match (a, b) {
            (None, None) => cmp::Ordering::Equal,
            (None, Some(_)) => cmp::Ordering::Greater,
            (Some(_), None) => cmp::Ordering::Less,
            (Some((_, aln_a, _)), Some((_, aln_b, _)))
                => aln_a.score.cmp(&aln_b.score)
        }
    })??;

    // previously, we've calculated this properly, using the alignment path
    // igblast instead just finds the start position, then goes 3 bases further
    // so we've mimicked their technique
    let cys_start_seq = find_index(&best_aln, &(cys_end_ref-3))?;
    let cys_end_seq = cys_start_seq + 3;

    // println!("best dist {} matched {} giving cys end seq {} which is ({}){}", 
    //     best_aln.score,
    //     name,
    //     cys_end_seq, 
    //     String::from_utf8(seq[cys_end_seq-3..cys_end_seq].to_vec()).expect("aa"),
    //     String::from_utf8(seq[cys_end_seq..cys_end_seq+10].to_vec()).expect("aa")
    // );

    if cys_end_seq >= seq.len() {
        // some reads get cut off just after the V gene (that is, FR3-IMGT)
        None
    } else {
        let (first_fr4_start, _, _) = fr4.clone()
            .find_all(&seq[cys_end_seq..], edit_dist)
            .min_by_key(|&(_, _, dist)| dist)?;
        let fr4_start = first_fr4_start + cys_end_seq;
        
        Some(&seq[cys_end_seq..fr4_start])
    }
}

fn find_index_old(ops: &[AlignmentOperation], ref_index: usize) -> Option<usize> {
    let mut seq_index = None;
    let (mut seq_i, mut ref_i) = (0, 0);

    for op in ops {
        if ref_i == ref_index {
            seq_index = Some(seq_i);
            // break
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

fn find_index(aln: &Alignment, ref_index: &usize) -> Option<usize> {
    let a = aln.path().into_iter().filter(|(x, _, _)| x.eq(ref_index));
    Some(a.last()?.1)
}

