# 🔎 find-cdr3

A fast program for extracting CDR3 regions from immunoglobulin sequences. [igblast](https://www.ncbi.nlm.nih.gov/igblast/) also does this, but it's slow and performs operations we're not interested in. This program is a stripped-down, sped-up version of igblast for our specific use case.

## Usage

Provide `find-cdr3` with an input fasta, a reference fasta, and an output destination. For each read in the input fasta, the CDR3 is located, and then the read's ID, sequence, and CDR3 sequence are printed to a line of an output CSV.

```
Usage: find-cdr3 [OPTIONS] --input-fasta <INPUT_FASTA> --reference-fasta <REFERENCE_FASTA> --output-csv <OUTPUT_CSV>

Options:
  -i, --input-fasta <INPUT_FASTA>
          Input fasta file of immunoglobin sequences
  -r, --reference-fasta <REFERENCE_FASTA>
          Reference fasta of different V-gene sequences. The Cys codon near the end of the V-gene sequence marks the start of the CDR3
  -o, --output-csv <OUTPUT_CSV>
          Output CSV
  -s, --sample-size <SAMPLE_SIZE>
          Size of pre-processing sample, pre-processed to find the most common reference sequences [default: 10000]
  -p, --parallel-chunk-size <PARALLEL_CHUNK_SIZE>
          Size of each chunk of reads to process in parallel [default: 500000]
  -e, --edit-dist <EDIT_DIST>
          Edit distance used for reference sequences [default: 10]
  -f, --fr4 <FR4>
          FR4 region. Located just after the CDR3 [default: TGGGGCAAAGGGACCCAGGTCAC]
  -h, --help
          Print help
```