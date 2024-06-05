#!/bin/bash +x

# Simulate the reads in fastq format.

../software/wgsim/wgsim -N 1000 -1 150 -2 150 -e 0.02 -r 0.05 -R 0.15 -X 0.3 ../src/hg19.fa.masked test_R1.fastq test_R2.fastq
