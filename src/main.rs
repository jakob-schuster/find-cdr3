use core::{panic, fmt};
use std::{io::{BufReader, BufRead, Write}, fs::File, cmp, string};

use bio::{pattern_matching::myers::{Myers, MyersBuilder}, io::fasta::Record, alphabets::protein, alignment::{Alignment, AlignmentOperation}};
use clap::{Parser, Arg};


#[derive(Parser,Debug)]
struct Args {
    #[arg(short,long)]
    input_fasta: String,
    #[arg(short,long)]
    reference_fasta: String,
    #[arg(short,long)]
    output_csv: String,
}

struct RefV {
    seq: String,
    cys_index: usize
}

#[cfg(test)]
mod tests {
    use super::*;
    fn test_translate() {
        assert_eq!(translate("TGGGGCCAGGGGACCCAGGTCACCGTC"), "WGQGTQVTV")
    }
}

fn translate(seq: &str) -> String {
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

fn last_cys(seq: &str) -> Option<usize> {
    let proteins = translate(seq);
    proteins.rfind('C').map_or(None, |a| Some(a*3))
}

fn parse_reference_b(reference_fasta: String) -> Vec<RefV> {
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

struct OutputRecord {
    name: String,
    seq: String,
    cdr3_seq: String
}

impl fmt::Display for OutputRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.name, self.seq, self.cdr3_seq)
    }
}

fn find_cdr3(seq: &String, reference_seqs: &Vec<RefV>, edit_dist: u8) -> Option<Vec<u8>> {
    // build ALL the myers matchers inside each iteration - bad?
    let reference_myers: Vec<(Myers, &usize)> = reference_seqs.into_iter()
        .map(|RefV { seq, cys_index } | (Myers::<u64>::new(seq.as_str().as_bytes()), cys_index))
        .collect();

    // first map to all the variable regions. get the best match (compare by edit dist)
    let v_matches = reference_myers.into_iter().map(|(mut myers, cys_index)| {
        let mut aln = Alignment::default();
        
        let mut matches = myers.find_all_lazy(seq.as_bytes(), edit_dist);
        let best_match = matches.by_ref().max_by_key(|&(_, dist)| dist);

        match best_match {
            None => None,
            Some((best_end, best_dist)) => {
                matches.alignment_at(best_end, &mut aln);
                
                // the start of the cys codon in the actual sequence, 
                // according to the position of cys in the reference 
                // and the operations involved in the alignment
                let cys_start_seq = aln.operations.into_iter()
                    .fold((0, 0), |(seq_i, v_i), operation| {
                        if v_i == *cys_index {
                            // already reached the cys; just pass it through unchanged
                            (seq_i, v_i)
                        } else {
                            // handle the current operation
                            match operation {
                                AlignmentOperation::Del => (seq_i, v_i + 1),
                                AlignmentOperation::Ins => (seq_i + 1, v_i),
                                AlignmentOperation::Match => (seq_i + 1, v_i + 1),
                                AlignmentOperation::Subst => (seq_i + 1, v_i + 1),
                                AlignmentOperation::Xclip(n) => (seq_i + n, v_i),
                                AlignmentOperation::Yclip(n) => (seq_i, v_i + n)
                            }
                        }
                    }
                ).0;
                
                Some((cys_start_seq, best_dist))
            }
        }
    });


    let v_best_match = v_matches.max_by(|a, b| {
        match (a, b) {
            (None, None) => cmp::Ordering::Equal,
            (None, Some(_)) => cmp::Ordering::Less,
            (Some(_), None) => cmp::Ordering::Greater,
            (Some((_, edit_dist_a)), Some((_, edit_dist_b))) => edit_dist_a.cmp(edit_dist_b) 
        }
    }).map_or(None, |a| {
        a.map_or(None, |(cys_start, _)| {
            Some(cys_start)
        })
    });
    
    
    // fold(None, |prev, m| {
    //     cmp::max(prev, m)
    // });

    // currently hardcoded codon for the fr4
    let mut fr4_imgt_myers = Myers::<u64>::new(b"TGGGGCAAAGGGACCCAGGTCAC");

    match v_best_match {
        None => None,
        Some(cys_start_seq) => {
            match fr4_imgt_myers.find_all(&seq.as_bytes()[cys_start_seq..], 3).min_by_key(|&(start, _, _)| start) {
                None => None,
                Some(first_fr4_match) => {
                    let fr4_start = first_fr4_match.0 + cys_start_seq;
                    Some(seq.as_bytes()[(cys_start_seq+3)..fr4_start].to_owned())
                }
            }
        } 
    }
}

fn parse_input(input_fasta: String, reference_seqs: &Vec<RefV>, edit_dist: u8) -> Vec<OutputRecord> {
    let reader = bio::io::fasta::Reader::new(
        BufReader::new(
            File::open(input_fasta).expect("Couldn't open input!")));

    let recs: Vec<_> = reader.records().collect();
    
    recs.into_iter()
        .map(|result| match result {
            Err(_) => panic!("Bad record in input!"),
            Ok(record) => {
                OutputRecord { 
                    name: record.id().to_string(),
                    seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                    cdr3_seq: match std::cmp::max(
                        find_cdr3(&String::from_utf8(record.seq().to_owned()).unwrap(), &reference_seqs, edit_dist), 
                        find_cdr3(&String::from_utf8(bio::alphabets::dna::revcomp(record.seq().to_owned())).unwrap(), &reference_seqs, edit_dist)) {

                        None => String::new(),
                        Some(bytes) => {String::from_utf8(bytes.to_vec()).unwrap()}
                    }
                }
            }
        }
    ).collect()
}

fn print_output(output_csv: String, output_records: &Vec<OutputRecord>) {
    let mut output_csv_file = File::create(output_csv)
        .expect("Couldn't create output file!");

    write!(output_csv_file, "id,sequence,cdr3_sequence\n")
        .expect("Couldn't write header line to output!");
    for record in output_records {
        write!(output_csv_file, "{}\n", record)
            .expect("Couldn't write line to output!");
    }
}

fn main() {
    let args = Args::parse();

    let reference_seqs_b = parse_reference_b(args.reference_fasta);
    let output_records = parse_input(args.input_fasta, &reference_seqs_b, 8);

    print_output(args.output_csv, &output_records);

    println!("hello")
}
