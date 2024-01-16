use std::{cmp, thread};

use bio::{alignment::{Alignment, AlignmentOperation}, pattern_matching::myers::Myers, io::fasta::Record};
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IndexedParallelIterator};

use crate::reference;

use super::OutputRecord;

pub(crate) fn parse_one_input(
    record: Record,
    reference_seqs: &Vec<reference::RefV>,
    edit_dist: u8,
    fr4: &Myers::<u64>
) -> OutputRecord {
    // println!("checking seq {}", record.id());

    let mut vec = Vec::new();
    
    vec![record.seq(), &bio::alphabets::dna::revcomp(record.seq())]
        .into_par_iter()
        .map(|seq| {
            find_cdr3(
                seq, 
                &reference_seqs.to_owned(),
                edit_dist,
                fr4
            )
        })
        .collect_into_vec(&mut vec);

    if let [forward, reverse] = &vec[..] {
        match (&forward, &reverse) {
            // neither was successful
            (None, None) => OutputRecord {
                name: String::from(record.id()),
                seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                cdr3_seq: String::from("neither")
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
                cdr3_seq: String::from("both")
            }
        }
    } else {
        panic!("Thread problem")
    }
}

fn find_index(ops: &[AlignmentOperation], ref_index: usize) -> Option<usize> {
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

fn find_index_new(aln: &Alignment, ref_index: &usize) -> Option<usize> {
    let mut a = aln.path().into_iter().filter(|(x, _, _)| x.eq(ref_index));
    let p = aln.path();

    Some(a.last()?.1)
}

pub(crate) fn find_cdr3(
    seq: &[u8], reference_seqs: &Vec<reference::RefV>, edit_dist: u8, fr4: &Myers::<u64>
) -> Option<Vec<u8>> {
    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_seqs.into_par_iter()
    .map(|reference::RefV { name, seq: _ref_seq, myers, cys_index: cys_end_ref } | {
        let mut aln = Alignment::default();
    
        let mut owned_myers = myers.clone();

        let mut matches = 
            owned_myers.find_all_lazy(seq, edit_dist);
        let (best_end, best_dist) = 
            matches.by_ref().min_by_key(|&(_, dist)| dist)?;

        matches.alignment_at(best_end, &mut aln);

        // let cys_ref_end = *cys_index + 3;
        // let cys_seq_end = find_index(&aln.clone().operations, cys_ref_end)?;
        
        Some((name, aln, cys_end_ref))
    });

    // get the cys start of the best match, if there is a match at all
    let (name, best_aln, cys_end_ref) = v_matches.min_by(|a, b| {
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
    let cys_start_seq = find_index_new(&best_aln, &(cys_end_ref-3))?;
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
        
        Some(seq[cys_end_seq..fr4_start].to_owned())
    }
}
