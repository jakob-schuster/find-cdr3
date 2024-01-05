# 🔎 find-cdr3

A fast program for extracting CDR3 regions from immunoglobulin sequences. [igblast](https://www.ncbi.nlm.nih.gov/igblast/) also does this, but it's significantly slower, and also performs lots of operations that we're not interested in. This program is intended as a stripped-down, sped-up version of igblast.

# Usage

Provide `find-cdr3` with an input fasta, a reference fasta, and an output destination.

```
Usage: find-cdr3 [OPTIONS] --input-fasta <INPUT_FASTA> --reference-fasta <REFERENCE_FASTA> --output-csv <OUTPUT_CSV>

Options:
  -i, --input-fasta <INPUT_FASTA>                  
  -r, --reference-fasta <REFERENCE_FASTA>          
  -o, --output-csv <OUTPUT_CSV>                    
  -s, --sample-size <SAMPLE_SIZE>                  [default: 10000]
  -p, --parallel-chunk-size <PARALLEL_CHUNK_SIZE>  [default: 500000]
  -n, --nondeterministic                           
  -e, --edit-dist <EDIT_DIST>                      [default: 10]
  -h, --help                                       Print help
```