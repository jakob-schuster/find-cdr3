use crate::reference::RefV;

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
                    // println!("new best is {} {} cys {} which is {}", best_edit_dist, ref_v.clone().seq, ref_v.clone().cys_index, &ref_v.seq[ref_v.cys_index..]);
                    // if it's the first, always update
                    Some((ref_v.clone(), best_edit_dist))
                },
                Some((_, old_edit_dist)) => {
                    if best_edit_dist < old_edit_dist {
                    // println!("new best is {} {} cys {} which is {}", best_edit_dist, ref_v.clone().name, ref_v.clone().cys_index, &ref_v.seq[ref_v.cys_index..]);
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