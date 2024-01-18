use bio;
use bio::pattern_matching::myers::{Myers, MyersBuilder};
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::find_v;
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct RefV {
    pub name: String,
    pub seq: String,
    pub myers: Myers,
    pub cys_index: usize,
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    #[test]
    fn test_translate() {
        assert_eq!(translate("TGGGGCCAGGGGACCCAGGTCACCGTC"), "WGQGTQVTV")
    }
}

pub(crate) fn translate(seq: &str) -> String {
    fn translate_codon(codon: &str) -> char {
        match codon {
            "TTT" => 'F',
            "TTC" => 'F',
            "TTA" => 'L',
            "TTG" => 'L',
            "TCT" => 'S',
            "TCC" => 'S',
            "TCA" => 'S',
            "TCG" => 'S',
            "TAT" => 'Y',
            "TAC" => 'Y',
            "TAA" => '-',
            "TAG" => '-',
            "TGT" => 'C',
            "TGC" => 'C',
            "TGA" => '-',
            "TGG" => 'W',

            "CTT" => 'L',
            "CTC" => 'L',
            "CTA" => 'L',
            "CTG" => 'L',
            "CCT" => 'P',
            "CCC" => 'P',
            "CCA" => 'P',
            "CCG" => 'P',
            "CAT" => 'H',
            "CAC" => 'H',
            "CAA" => 'Q',
            "CAG" => 'Q',
            "CGT" => 'R',
            "CGC" => 'R',
            "CGA" => 'R',
            "CGG" => 'R',

            "ATT" => 'I',
            "ATC" => 'I',
            "ATA" => 'I',
            "ATG" => 'M',
            "ACT" => 'T',
            "ACC" => 'T',
            "ACA" => 'T',
            "ACG" => 'T',
            "AAT" => 'N',
            "AAC" => 'N',
            "AAA" => 'K',
            "AAG" => 'K',
            "AGT" => 'S',
            "AGC" => 'S',
            "AGA" => 'R',
            "AGG" => 'R',

            "GTT" => 'V',
            "GTC" => 'V',
            "GTA" => 'V',
            "GTG" => 'V',
            "GCT" => 'A',
            "GCC" => 'A',
            "GCA" => 'A',
            "GCG" => 'A',
            "GAT" => 'D',
            "GAC" => 'D',
            "GAA" => 'E',
            "GAG" => 'E',
            "GGT" => 'G',
            "GGC" => 'G',
            "GGA" => 'G',
            "GGG" => 'G',

            _ => '?',
        }
    }

    if seq.len() < 3 {
        // terminate on the end of sequences
        String::new()
    } else {
        // or, recursively call
        format!("{}{}", translate_codon(&seq[..3]), translate(&seq[3..]))
    }
}

pub(crate) fn last_cys(seq: &str) -> Option<usize> {
    let proteins = translate(seq);
    proteins.rfind('C').map_or(None, |a| Some(a * 3))
}

pub(crate) fn last_cys_simple(seq: &str) -> Option<usize> {
    let mut cys_myers = Myers::<u64>::new(b"TGT");
    // let mut cys_myers = MyersBuilder::new().ambig(b'B', b"TC")
    //     .build_64(b"TGB");

    cys_myers
        .find_all_end(seq.as_bytes(), 0)
        .max_by_key(|(st, ed)| st.clone())
        .map_or(None, |(st, ed)| Some(st + 1))
}

pub(crate) fn parse_reference(reference_fasta: &str) -> Vec<RefV> {
    let reader = bio::io::fasta::Reader::new(BufReader::new(
        File::open(reference_fasta).expect("Couldn't open reference!"),
    ));

    reader
        .records()
        .into_iter()
        .map(|result| match result {
            Ok(record) => {
                let seq = String::from_utf8(record.seq().to_vec())
                    .expect("Bad sequence line!")
                    .to_ascii_uppercase();

                let short_seq_start = cmp::max(0, seq.len() - 64);
                let short_seq = &seq[short_seq_start..];

                RefV {
                    name: String::from(record.id()),
                    seq: String::from(short_seq),
                    myers: Myers::<u64>::new(short_seq.as_bytes()),
                    cys_index: last_cys(&seq).expect("No Cys!") - short_seq_start + 3,
                }
            }
            Err(_) => panic!("Bad record!"),
        })
        .collect()
}

pub fn optimise_refs(
    reference_seqs: &Vec<RefV>,
    records: bio::io::fasta::Records<BufReader<Box<dyn BufRead>>>,
    sample_size: usize,
    parallel_chunk_size: usize,
    edit_dist: u8,
    reference_size: usize,
) -> (Vec<RefV>, Vec<RefV>) {
    // go through the records populating a counts map
    let mut map: HashMap<RefV, usize> = HashMap::new();

    for bunch in &records.take(sample_size).chunks(parallel_chunk_size) {
        let mut outs = Vec::new();

        bunch
            .collect_vec()
            .into_par_iter()
            .map(|result| match result {
                Err(_) => panic!("Bad record!"),
                Ok(record) => find_v::parse_one(record.seq(), &reference_seqs, true, edit_dist),
            })
            .collect_into_vec(&mut outs);

        for out in outs {
            match out {
                None => {}
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

    let sorted_map = map
        .into_iter()
        .sorted_by(|(_, edit_a), (_, edit_b)| edit_a.cmp(edit_b))
        .rev()
        .collect_vec();

    for (k, v) in &sorted_map {
        println!("{},{}", k.name, v);
    }

    // convert the map into a sorted list
    let mut best = sorted_map
        .into_iter()
        .map(|(ref_v, _edit_dist)| ref_v)
        .collect_vec();

    let mut rest = reference_seqs
        .into_iter()
        .filter(|&a| -> bool { best.contains(a) })
        .map(|a| a.to_owned())
        .collect_vec();

    best.append(&mut rest);

    // best
    (
        best.clone().into_iter().take(reference_size).collect_vec(),
        best.into_iter().dropping(reference_size).collect_vec(),
    )
}
