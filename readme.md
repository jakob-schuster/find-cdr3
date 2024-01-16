# 🔎 find-cdr3

A fast program for extracting CDR3 regions from immunoglobulin sequences. [igblast](https://www.ncbi.nlm.nih.gov/igblast/) also does this, but it's slow and performs operations we're not interested in. This program is a stripped-down, sped-up version of igblast for our specific use case. Around 91.4% of `find-cdr3`'s results perfectly match igblast; the rest are only very slight variations due to different 

## Usage

Provide `find-cdr3` with an input fasta, a reference fasta, and an output destination. For each read in the input fasta, the CDR3 is located, and then the read's ID, sequence, and CDR3 sequence are printed to a line of an output TSV.

```
Usage: find-cdr3 [OPTIONS] --reference-fasta <REFERENCE_FASTA>

Options:
  -i, --input-fasta <INPUT_FASTA>
          Input fasta file of immunoglobin sequences. Reads from stdin if none is provided [default: stdin]
  -r, --reference-fasta <REFERENCE_FASTA>
          Reference fasta of different V-gene sequences. The Cys codon near the end of the V-gene sequence marks the start of the CDR3
  -o, --output-tsv <OUTPUT_TSV>
          Output TSV file of reads and their CDR3s. Writes to stdout if none is provided [default: stdout]
      --sample-size <SAMPLE_SIZE>
          Size of pre-processing sample, pre-processed to find the most common reference sequences. Not used currently [default: 10000]
      --reference-size <REFERENCE_SIZE>
          Number of top V-gene sequences to actually search for. The less you use, the quicker the program, but also the less accurate [default: 20]
  -c, --parallel-chunk-size <PARALLEL_CHUNK_SIZE>
          Size of each chunk of reads to process in parallel [default: 500000]
  -e, --edit-dist <EDIT_DIST>
          Edit distance used for reference sequences [default: 20]
  -f, --fr4 <FR4>
          FR4 region. Located just after the CDR3 [default: TGGGGCAAAGGGACCCAGGTCAC]
  -l, --headerless
          Omit the header from the output. Useful if you're piping into another program to do more processing
  -h, --help
          Print help
```