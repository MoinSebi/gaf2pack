# gaf2pack

Convert GAF alignments (https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) to pack format. 
Processing cigar string for additional information. Multithreading via rayon. Pack (coverage) output is base-level. Run it with multiple threads to speed up the process otherwise it will take a long time.


## Installation

```
git clone https://github.com/MoinSebi/gaf2pack
cd gaf2pack
cargo build --release
./target/release/gaf2pack -h 
```

## Usage 

```
./gaf2pack -g input.gfa -o output.pack -a alignment.gaf -t 10 
```