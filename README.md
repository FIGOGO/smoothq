## Introduction
Smoothq is a sensitive overlap detection program 
that detects overlaps for long but error-prone DNA reads
generated by the third generation sequencing technology 
such as PacBio SMRT and Oxford Nanopore.

The innovation of Smoothq is to convert q-gram to smooth q-gram via 
[CGK-embedding][CGK] and resampling,
which can capture q-gram paris with small edit distance.
Using smooth q-gram instead of q-gram as seeds, 
Smoothq is very sensitive to detect overlaps 
between DNA reads even when their overlap length is small.

## Getting Started
### Compile
```sh
git clone git@github.com:FIGOGO/smoothq.git
cd smoothq
cmake .
make
```

### Usage
General usage:
```sh
./smoothq input.fasta [options] > overlap.txt
```

A small dataset *ecoli1000.fasta* is provided for testing purposes:

```sh
./smoothq ecoli1000.fasta -t 4 -q 16 > overlap.txt
```

## Options
The size of the q-gram, cgk-embedding and smooth q-gram can be tune by the following:
    
- -q, default = 14: q-gram size
- -e, default = 35: embedding size
- -m, default = 16: smooth q-gram size

Parameters controlling the number of signatures (sensitivity) can be tuned by the following:

- -f, default = 5e-4: filtering threshold
- -c, default = 5:    minimum number of signature match
- -r, default = 0.2:  sampling rate

Number of threads used in the algorithm can be tuned by the following:
- -t, default = 16: number of threads

## Output format
Our output is similar to BLASR’s M4 format.

1. Index of sequence 1 (Index is 1-based)
2. Index of sequence 2
3. Sequence tags in original fasta file separated by '_'
4. Signature density of the overlap
5. Direction of sequence 1 (0 = forward, 1 = reversed)
6. Overlap start position on sequence 1
7. Overlap end position on sequence 1
8. Length of sequence 1
9. Direction of sequence 2 (0 = forward, 1 = reversed)
10. Overlap start position on sequence 2
11. Overlap end position on sequence 2
12. Length of sequence 2

## Acknowledgement
We use the [xxHash][XHS] package in smooth q-gram hashing. 

[XHS]: https://github.com/Cyan4973/xxHash
[CGK]: https://dl.acm.org/citation.cfm?id=2897577
