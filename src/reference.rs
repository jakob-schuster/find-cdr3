use std::cmp;
use std::fs::File;
use std::io::BufReader;
use bio;

pub(crate) struct RefV {
    pub(crate) seq: String,
    pub(crate) cys_index: usize
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    fn test_translate() {
        assert_eq!(translate("TGGGGCCAGGGGACCCAGGTCACCGTC"), "WGQGTQVTV")
    }
}

pub(crate) fn translate(seq: &str) -> String {
    match seq {
        "TTT" => "F".to_string(), "TTC" => "F".to_string(), "TTA" => "L".to_string(), "TTG" => "L".to_string(),
        "TCT" => "S".to_string(), "TCC" => "S".to_string(), "TCA" => "S".to_string(), "TCG" => "S".to_string(),
        "TAT" => "Y".to_string(), "TAC" => "Y".to_string(), "TAA" => "-".to_string(), "TAG" => "-".to_string(),
        "TGT" => "C".to_string(), "TGC" => "C".to_string(), "TGA" => "-".to_string(), "TGG" => "W".to_string(),

        "CTT" => "L".to_string(), "CTC" => "L".to_string(), "CTA" => "L".to_string(), "CTG" => "L".to_string(),
        "CCT" => "P".to_string(), "CCC" => "P".to_string(), "CCA" => "P".to_string(), "CCG" => "P".to_string(),
        "CAT" => "H".to_string(), "CAC" => "H".to_string(), "CAA" => "Q".to_string(), "CAG" => "Q".to_string(),
        "CGT" => "R".to_string(), "CGC" => "R".to_string(), "CGA" => "R".to_string(), "CGG" => "R".to_string(),

        "ATT" => "I".to_string(), "ATC" => "I".to_string(), "ATA" => "I".to_string(), "ATG" => "M".to_string(),
        "ACT" => "T".to_string(), "ACC" => "T".to_string(), "ACA" => "T".to_string(), "ACG" => "T".to_string(),
        "AAT" => "N".to_string(), "AAC" => "N".to_string(), "AAA" => "K".to_string(), "AAG" => "K".to_string(),
        "AGT" => "S".to_string(), "AGC" => "S".to_string(), "AGA" => "R".to_string(), "AGG" => "R".to_string(),

        "GTT" => "V".to_string(), "GTC" => "V".to_string(), "GTA" => "V".to_string(), "GTG" => "V".to_string(),
        "GCT" => "A".to_string(), "GCC" => "A".to_string(), "GCA" => "A".to_string(), "GCG" => "A".to_string(),
        "GAT" => "D".to_string(), "GAC" => "D".to_string(), "GAA" => "E".to_string(), "GAG" => "E".to_string(),
        "GGT" => "G".to_string(), "GGC" => "G".to_string(), "GGA" => "G".to_string(), "GGG" => "G".to_string(),

        _ => if seq.len() < 3 {
            "?".to_string()
        } else {
            format!("{}{}", translate(&seq[0..3]), translate(&seq[3..]))
        }
    }
}

pub(crate) fn last_cys(seq: &str) -> Option<usize> {
    let proteins = translate(seq);
    proteins.rfind('C').map_or(None, |a| Some(a*3))
}

pub(crate) fn parse_reference(reference_fasta: String) -> Vec<RefV> {
    let reader = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(reference_fasta).expect("Couldn't open reference!")));

    reader.records().into_iter().map(|result| {
        match result {
            Ok(record) => {String::from_utf8(record.seq()[cmp::max(0,record.seq().len()-63)..].to_ascii_uppercase().to_vec()).expect("Bad sequence line!")}
            Err(_) => panic!("Bad record!")
        }
    }).map(|seq| {RefV { seq: seq.clone(), cys_index: last_cys(seq.as_str()).expect("No Cys") } }).collect()
}