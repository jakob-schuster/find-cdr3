# 🔎 find-cdr3

A fast program for extracting CDR3 regions from immunoglobulin sequences. [igblast](https://www.ncbi.nlm.nih.gov/igblast/) also does this, but it's slow and performs operations we're not interested in. This program is a stripped-down, hopefully sped-up version of igblast for our specific use case. Around 91.4% of `find-cdr3`'s results perfectly match igblast's; the rest are only very slight variations due to different alignment path strategies.

Download latest nightly builds for [Linux](https://nightly.link/jakob-schuster/find-cdr3/workflows/rust/main/find-cdr3-x86_64-unknown-linux-gnu.zip) / [Mac](https://nightly.link/jakob-schuster/find-cdr3/workflows/rust/main/find-cdr3-x86_64-apple-darwin.zip) / [Windows](https://nightly.link/jakob-schuster/find-cdr3/workflows/rust/main/find-cdr3-x86_64-pc-windows-msvc.zip)

## 🧭 Usage

Provide `find-cdr3` with an input fasta, a reference fasta, and an output destination. For each read in the input fasta, the CDR3 is located, and then the read's ID, sequence, and CDR3 sequence are printed to a line of an output TSV.

```
Usage: find-cdr3 [OPTIONS] --reference-fasta <REFERENCE_FASTA>

Options:
  -i, --input-reads <INPUT_READS>
          Input fasta/fastq file of immunoglobin sequences. Reads from stdin if none is provided [default: stdin]
  -r, --reference-fasta <REFERENCE_FASTA>
          Reference fasta of different V-gene sequences. The Cys codon near the end of the V-gene sequence marks the start of the CDR3
  -o, --output-tsv <OUTPUT_TSV>
          Output TSV file of reads and their CDR3s. Writes to stdout if none is provided [default: stdout]
      --sample-size <SAMPLE_SIZE>
          Size of pre-processing sample, pre-processed to find the most common reference sequences [default: 10000]
      --reference-size <REFERENCE_SIZE>
          Number of top V-gene sequences to actually search for. The less you use, the quicker the program, but also the less accurate [default: 5]
  -t, --threads <THREADS>
          Number of threads to use for parallel processing [default: 1]
  -c, --parallel-chunk-size <PARALLEL_CHUNK_SIZE>
          Size of each chunk of reads to process in parallel [default: 300000]
  -e, --edit-dist <EDIT_DIST>
          Edit distance used for reference sequences [default: 20]
  -f, --fr4 <FR4>
          FR4 region. Located just after the CDR3 [default: TGGGGCAAAGGGACCCAGGTCAC]
  -l, --headerless
          Omit the header from the output. Useful if you're piping into another program to do more processing
  -h, --help
          Print help
```

## 🔩 Implementation

The `rust-bio` crate was used for parsing and matching patterns in reads. All alignment is done using Myers bit-vector algorithm.

### Preprocessing

The reference file (`reference-fasta`) is read first. For each reference sequence, the Cys codon's position is calculated.

Since there are often hundreds of reference V-genes of which only a few will be present in a given dataset, we do some preprocessing to cull the less relevant references. The first `sample-size` reads are taken from `input-fasta`, and for each read, the best matching V-gene is recorded. Then, the V-genes are sorted and only the `reference-size` most frequent references are kept for processing. Depending on your data, this may be a lossy optimisation - to avoid it, increase `reference-size` to be large enough for all of your references to be kept.

### Processing

The input file `input-reads` is parsed. `parallel-chunk-size` reads are processed, using `threads` threads, and then the results are written to an output file by one thread, before the next chunk is handled. 

Each read is scanned for a V-gene match allowing for `edit-dist` mismatches, selecting the one with the fewest mismatches. Then, based on the position of the Cys codon in the reference sequence and the path of alignment between the reference sequence and the read, the location of the Cys codon within the read is calculated. Then, the FR4 (`fr4`) region is located, allowing for `edit-dist` mismatches. If both the V-gene and the FR4 region could be found, the region from the Cys codon to the start of the FR4 is taken as the CDR3, and extracted. 

For each read, the TSV output file contains: 
- The read ID
- The full sequence of the read 
- The name of the V-gene 
- The CDR3 sequence 
- The translated amino acid sequence spanning from the start of the V-gene to the end of the FR4

## 🗺 Todo

- Add file splitting and merging for further speed improvement