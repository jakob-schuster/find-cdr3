use std::cmp;

use bio::{pattern_matching::myers::Myers, alignment::Alignment};
use input::Read;
use itertools::Itertools;
use util::Ran;

use crate::reference;
use super::*;

/// A single result of the `parse_read` computation.
/// Contains all the information we want to extract,
/// and all the types of read we categorise.
#[derive(Default, Debug)]
pub enum OutputRecord {
    #[default]
    Neither,
    Forward(String, Vec<u8>, Vec<u8>),
    Reverse(String, Vec<u8>, Vec<u8>),
    Both
}

impl OutputRecord {
    pub fn fmt_seq(&self) -> String {
        match self {
            OutputRecord::Neither => 
                String::from("\t\t\t"),
            OutputRecord::Forward(name, cdr3, full_seq) =>
                format!("{}\t{}\t{}",
                    name,
                    std::str::from_utf8(cdr3).unwrap(),
                    crate::reference::translate(full_seq),
                ),
            OutputRecord::Reverse(name, cdr3, full_seq) =>
                format!("{}\t{}\t{}",
                    name,
                    std::str::from_utf8(cdr3).unwrap(),
                    crate::reference::translate(full_seq),
                    // crate::reference::translate(&seq[full_seq.start..full_seq.end])
                ),
            OutputRecord::Both =>
                String::from("\t\t\t"),
        }
    }
}

/// Parses a single read, trying to find and extract the cdr3 region 
/// and all other information. Tries in both read orientations.
pub fn parse_read(
    read: &Read,
    reference_seqs: &[reference::RefV],
    edit_dist: u8,
    fr4: &Myers::<u64>
) -> OutputRecord {
    let seq = read.seq();
    let rev_seq = &bio::alphabets::dna::revcomp(seq);

    let vec = [seq, rev_seq]
        .iter()
        .map(|seq| {
            find_cdr3(seq, reference_seqs, edit_dist, fr4)
        }).collect_vec();

    match &vec[..] {
        // neither was successful
        [None, None] => OutputRecord::Neither,

        // forwards was successful
        [Some((name, seq_cdr3, full_seq)), None] => 
            OutputRecord::Forward(name.clone(), seq_cdr3.clone(), full_seq.clone()),

        // reverse was successful
        [None, Some((name, rev_seq_cdr3, rev_full_seq))] => 
            OutputRecord::Reverse(
                name.clone(), rev_seq_cdr3.clone(), rev_full_seq.clone()
            ),

        // both were successful - ambiguous
        [Some(_), Some(_)] => OutputRecord::Both,

        _ => panic!("Thread problem")
    }
}

/// Searches through a single sequence, trying to find and extract the cdr3
/// region. 
fn find_cdr3(
    seq: &[u8], 
    reference_seqs: &[reference::RefV], 
    edit_dist: u8, 
    fr4: &Myers::<u64>
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
    let cys_start_seq = find_cys_start_index(&best_aln, &(cys_end_ref-3))?;
    let cys_end_seq = cys_start_seq + 3;

    if cys_end_seq >= seq.len() {
        // some reads get cut off just after the V gene (that is, FR3-IMGT)
        None
    } else {
        let (best_fr4_start, best_fr4_end, _) = fr4.clone()
            .find_all(&seq[cys_end_seq..], edit_dist)
            .min_by_key(|&(_, _, dist)| dist)?;
        let fr4_start = best_fr4_start + cys_end_seq;
        
        Some((
            v_name.clone(), 
            Vec::from(&seq[cys_end_seq..fr4_start]), 
            Vec::from(&seq[best_aln.ystart..best_fr4_end + cys_end_seq])
        ))
    }
}

/// Taking an alignment path for the V gene, finds the index of the Cys start within the sequence.
/// This is considered the start of the cdr3.
fn find_cys_start_index(aln: &Alignment, ref_index: &usize) -> Option<usize> {
    let a = aln.path().into_iter().filter(|(x, _, _)| x.eq(ref_index));
    Some(a.last()?.1)
}

