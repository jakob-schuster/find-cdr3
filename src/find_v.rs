use std::{cmp, io::BufReader, fs::File, collections::{HashMap, HashSet}};

use bio::{io::fasta::{Record, Records}, pattern_matching::myers::Myers, alignment::Alignment};
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator, IndexedParallelIterator};

use crate::reference::{RefV, self};

/// Takes one record, aligns it to all the V-genes, 
/// returns the one against which it matched
pub(crate) fn parse_one(
    seq: &[u8],
    reference_seqs: &Vec<RefV>,
    check_reverse: bool,
    edit_dist: u8
) -> Option<(RefV, u8)> {
    let mut best_ref_v: Option<(RefV, u8)> = None;
    for ref_v in reference_seqs {
        let (_, best_edit_dist) = ref_v.myers.find_best_end(seq);

        if best_edit_dist < edit_dist {
            // it is below the threshold and worth considering
            best_ref_v = match best_ref_v {
                None => {
                    // if it's the first, always update
                    Some((ref_v.clone(), best_edit_dist))
                },
                Some((_, old_edit_dist)) => {
                    if best_edit_dist < old_edit_dist {
                        // if it's better, update
                        Some((ref_v.clone(), best_edit_dist))
                    } else {
                        // otherwise, leave unchanged
                        best_ref_v
                    }
                }
            }
        }
    }

    if check_reverse {
        // may need to do annoying pattern match check on the reverse one
        match (&best_ref_v, &parse_one(
            &bio::alphabets::dna::revcomp(seq), 
            reference_seqs, 
            false, 
            edit_dist
        )) {
            (a, None) => a.clone(),
            (None, b) => b.clone(),
            (a@Some((_, edit_dist_a)), b@Some((_, edit_dist_b))) => {
                // select the one with the lowest edit distance
                if edit_dist_a < edit_dist_b { a.clone() } else { b.clone() }
            }
        }
    } else {
        // otherwise, just go with this
        best_ref_v.clone()
    }
}

/// Takes one record, aligns it to all the V-genes, 
/// returns the one against which it matched
pub(crate) fn parse_one_lazy(
    seq: &[u8],
    reference_seqs: &Vec<RefV>,
    edit_dist: u8
) 
-> Option<(RefV, u8)> 
{
    let mut aln = Alignment::default();
    for ref_v in reference_seqs {
        let mut owned_myers = ref_v.myers.clone();

        let mut matches = 
            owned_myers.find_all_lazy(seq, edit_dist);
        let best_match = 
            matches.by_ref().max_by_key(|&(_, dist)| dist);

        match best_match {
            None => {}
            Some((best_end, best_dist)) => {
                
            }
        }
        if edit_dist > 3 {
            
            // it is below the threshold and worth considering
            // best_ref_v = match best_ref_v {
            //     None => {
            //         // if it's the first, always update
            //         Some((ref_v.clone(), best_edit_dist))
            //     },
            //     Some((_, old_edit_dist)) => {
            //         if best_edit_dist < old_edit_dist {
            //             // if it's better, update
            //             Some((ref_v.clone(), best_edit_dist))
            //         } else {
            //             // otherwise, leave unchanged
            //             best_ref_v
            //         }
            //     }
            // }
        }

    }
    return None
}

pub fn optimise_refs(reference_seqs: &Vec<RefV>, records: Records<BufReader<BufReader<File>>>, sample_size: usize, parallel_chunk_size: usize, edit_dist: u8) -> Vec<RefV> {
    // go through the records populating a counts map
    let mut map: HashMap<reference::RefV, usize> = HashMap::new();

    for bunch in &records.take(sample_size).chunks(parallel_chunk_size) {
        let mut outs = Vec::new();
        
        bunch.collect_vec().into_par_iter().map(|result| {
            match result {
                Err(_) => panic!("Bad record!"),
                Ok(record) => parse_one(record.seq(), &reference_seqs, true, edit_dist)
            }
        }).collect_into_vec(&mut outs);

        for out in outs {
            match out {
                None => { } 
                Some((ref_v, _)) => {
                    let new_size = match map.get(&ref_v) {
                        Some(size) => size + 1,
                        None => 1,
                    };
                    map.insert(ref_v, new_size);
                }
            }
        }    
    }

    let sorted_map = map.into_iter()
        .sorted_by(|a@(_, edit_a), b@(_, edit_b)| {edit_a.cmp(edit_b)})
        .rev()
        .collect_vec();

    for (k, v) in &sorted_map {
        println!("{},{}", k.seq, v);
    }

    // convert the map into a sorted list
    let mut best = sorted_map.into_iter()
        .map(|(ref_v, edit_dist)| { ref_v })
        .collect_vec();

    let mut rest = reference_seqs.into_iter()
        .filter(|&a| -> bool {best.contains(a)})
        .map(|a| a.to_owned())
        .collect_vec();

    best.append(&mut rest);
    
    // best
    best.into_iter().take(15).collect_vec()
}