use std::cmp;
use std::fs::File;
use std::io::BufReader;
use bio;
use bio::pattern_matching::myers::Myers;
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct RefV {
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

pub(crate) fn parse_reference(reference_fasta: &str) -> Vec<RefV> {
    let reader = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(reference_fasta).expect("Couldn't open reference!")));

    reader.records().into_iter().map(|result| {
        match result {
            Ok(record) => {
                String::from_utf8(
                    record.seq()[cmp::max(0,record.seq().len()-63)..]
                    .to_ascii_uppercase()
                ).expect("Bad sequence line!")
            }
            Err(_) => panic!("Bad record!")
        }
    }).map(|seq| {RefV {
        seq: seq.to_owned(),
        myers: Myers::<u64>::new(seq.as_bytes()),
        cys_index: last_cys(&seq).expect("No Cys") 
    }}).collect()
}