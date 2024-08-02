use bio::pattern_matching::myers::Myers;
use itertools::Itertools;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use crate::{find_v, Args};
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct RefV {
    pub name: String,
    pub seq: Vec<u8>,
    pub myers: Myers,
    pub cys_index: usize,
}

/// Translates a sequence into protein.
pub fn translate(seq: &[u8]) -> String {
    fn translate_codon(codon: &[u8]) -> char {
        match codon {
            b"TTT" => 'F',
            b"TTC" => 'F',
            b"TTA" => 'L',
            b"TTG" => 'L',
            b"TCT" => 'S',
            b"TCC" => 'S',
            b"TCA" => 'S',
            b"TCG" => 'S',
            b"TAT" => 'Y',
            b"TAC" => 'Y',
            b"TAA" => '-',
            b"TAG" => '-',
            b"TGT" => 'C',
            b"TGC" => 'C',
            b"TGA" => '-',
            b"TGG" => 'W',

            b"CTT" => 'L',
            b"CTC" => 'L',
            b"CTA" => 'L',
            b"CTG" => 'L',
            b"CCT" => 'P',
            b"CCC" => 'P',
            b"CCA" => 'P',
            b"CCG" => 'P',
            b"CAT" => 'H',
            b"CAC" => 'H',
            b"CAA" => 'Q',
            b"CAG" => 'Q',
            b"CGT" => 'R',
            b"CGC" => 'R',
            b"CGA" => 'R',
            b"CGG" => 'R',

            b"ATT" => 'I',
            b"ATC" => 'I',
            b"ATA" => 'I',
            b"ATG" => 'M',
            b"ACT" => 'T',
            b"ACC" => 'T',
            b"ACA" => 'T',
            b"ACG" => 'T',
            b"AAT" => 'N',
            b"AAC" => 'N',
            b"AAA" => 'K',
            b"AAG" => 'K',
            b"AGT" => 'S',
            b"AGC" => 'S',
            b"AGA" => 'R',
            b"AGG" => 'R',

            b"GTT" => 'V',
            b"GTC" => 'V',
            b"GTA" => 'V',
            b"GTG" => 'V',
            b"GCT" => 'A',
            b"GCC" => 'A',
            b"GCA" => 'A',
            b"GCG" => 'A',
            b"GAT" => 'D',
            b"GAC" => 'D',
            b"GAA" => 'E',
            b"GAG" => 'E',
            b"GGT" => 'G',
            b"GGC" => 'G',
            b"GGA" => 'G',
            b"GGG" => 'G',

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

/// Finds the last Cys codon in a sequence.
fn last_cys(seq: &[u8]) -> Option<usize> {
    let proteins = translate(seq);
    proteins.rfind('C')
        .map(|a| a * 3)
}

/// Parses the reference. Each reference sequence is trucated to <=64 
/// characters to fit in a rust-bio Myers bit-vector (adjusted to 
/// be in the correct coding frame). The position of the last Cys
/// codon is also found and recorded.
pub fn parse_reference(reference_fasta: &str) -> Vec<RefV> {
    let reader = bio::io::fasta::Reader::new(BufReader::new(
        File::open(reference_fasta).expect("Couldn't open reference!"),
    ));

    reader
        .records()
        .map(|result| match result {
            Ok(record) => {
                let seq = record.seq()
                    .to_ascii_uppercase();

                let short_seq_len = 64;
                let reading_frame = (seq.len() - short_seq_len) % 3;
                let short_seq_start = cmp::max(0, seq.len() - short_seq_len + (3 - reading_frame));

                let short_seq = &seq[short_seq_start..];

                RefV {
                    name: String::from(record.id()),
                    seq: Vec::from(short_seq),
                    myers: Myers::<u64>::new(short_seq),
                    cys_index: last_cys(&seq).expect("No Cys!") - short_seq_start + 3,
                }
            }
            Err(_) => panic!("Bad record!"),
        })
        .collect()
}

/// Takes a sample of the input, finding the reference_size most common V genes
/// and dropping the rest.
pub fn optimise_refs(
    reference_seqs: &Vec<RefV>,
    args: &Args,
) -> (Vec<RefV>, Vec<RefV>) {
    // go through the records populating a counts map
    let mut map: HashMap<RefV, usize> = HashMap::new();

    let reader = crate::input::Reader::new(
        &args.input_reads, 
        args.parallel_chunk_size, 
        Some(args.sample_size)
    ).unwrap();
    
    reader.map(
        |read| {
            find_v::parse_one(read.seq(), reference_seqs, true, args.edit_dist)
        },
        |opt, map: &mut HashMap<RefV, usize>| {
            if let Some((ref_v, _)) = opt {
                let new_size = match map.get(ref_v) {
                    Some(size) => size + 1,
                    None => 1
                };

                map.insert(ref_v.clone(), new_size);
            }
        },
        &mut map
    );

    // sort the map to find the best V genes
    let sorted_map = map
        .iter()
        .sorted_by(|(_, edit_a), (_, edit_b)| edit_a.cmp(edit_b))
        .map(|(ref_v, _)| ref_v.clone())
        .rev();

    // only take the best reference_size results
    let best = sorted_map.clone()
        .take(args.reference_size);
    let rest = sorted_map.clone()
        .dropping(args.reference_size);
    
    (best.collect(), rest.collect())
}


#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    #[test]
    fn test_translate() {
        assert_eq!(translate(b"TGGGGCCAGGGGACCCAGGTCACCGTC"), "WGQGTQVTV")
    }
}