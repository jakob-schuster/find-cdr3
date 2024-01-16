use std::cmp;
use std::fs::File;
use std::io::BufReader;
use bio;
use bio::pattern_matching::myers::{Myers, MyersBuilder};
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct RefV {
    pub name: String,
    pub seq: String,
    pub myers: Myers,
    pub cys_index: usize
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    fn test_translate() {
        assert_eq!(translate("TGGGGCCAGGGGACCCAGGTCACCGTC"), "WGQGTQVTV")
    }
}

pub(crate) fn translate(seq: &str) -> String {
    fn translate_codon(codon: &str) -> char {
        match codon {
            "TTT" => 'F', "TTC" => 'F', "TTA" => 'L', "TTG" => 'L',
            "TCT" => 'S', "TCC" => 'S', "TCA" => 'S', "TCG" => 'S',
            "TAT" => 'Y', "TAC" => 'Y', "TAA" => '-', "TAG" => '-',
            "TGT" => 'C', "TGC" => 'C', "TGA" => '-', "TGG" => 'W',
    
            "CTT" => 'L', "CTC" => 'L', "CTA" => 'L', "CTG" => 'L',
            "CCT" => 'P', "CCC" => 'P', "CCA" => 'P', "CCG" => 'P',
            "CAT" => 'H', "CAC" => 'H', "CAA" => 'Q', "CAG" => 'Q',
            "CGT" => 'R', "CGC" => 'R', "CGA" => 'R', "CGG" => 'R',
    
            "ATT" => 'I', "ATC" => 'I', "ATA" => 'I', "ATG" => 'M',
            "ACT" => 'T', "ACC" => 'T', "ACA" => 'T', "ACG" => 'T',
            "AAT" => 'N', "AAC" => 'N', "AAA" => 'K', "AAG" => 'K',
            "AGT" => 'S', "AGC" => 'S', "AGA" => 'R', "AGG" => 'R',
    
            "GTT" => 'V', "GTC" => 'V', "GTA" => 'V', "GTG" => 'V',
            "GCT" => 'A', "GCC" => 'A', "GCA" => 'A', "GCG" => 'A',
            "GAT" => 'D', "GAC" => 'D', "GAA" => 'E', "GAG" => 'E',
            "GGT" => 'G', "GGC" => 'G', "GGA" => 'G', "GGG" => 'G',
    
            _ => '?'
        }
    }

    if seq.len() < 3 { // terminate on the end of sequences
        String::new()
    } else { // or, recursively call
        format!("{}{}", translate_codon(&seq[..3]), translate(&seq[3..]))
    }
}

pub(crate) fn last_cys(seq: &str) -> Option<usize> {
    let proteins = translate(seq);
    proteins.rfind('C').map_or(None, |a| Some(a*3))
}

pub(crate) fn last_cys_simple(seq: &str) -> Option<usize> {
    let mut cys_myers = Myers::<u64>::new(b"TGT");
    // let mut cys_myers = MyersBuilder::new().ambig(b'B', b"TC")
    //     .build_64(b"TGB");

    cys_myers.find_all_end(seq.as_bytes(), 0)
        .max_by_key(|(st, ed)| st.clone())
        .map_or(None, |(st, ed)| Some(st + 1))
}

pub(crate) fn parse_reference(reference_fasta: &str) -> Vec<RefV> {
    let reader = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(reference_fasta).expect("Couldn't open reference!")));

    reader.records().into_iter().map(|result| {        
        match result {
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
            Err(_) => panic!("Bad record!")
        }
    }).collect()
}