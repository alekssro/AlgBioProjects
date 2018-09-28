#!/bin/bash
# Script to measure the running time of our implementation of global_alignment.py

test_seq1="tatggagagaataaaagaactgagagatctaatgtcgcagtcccgcactcgcgagatact"
test_seq2="tccaaaatggaagactttgtgcgacaatgcttcaatccaatgatcgtcgagcttgcggaa"
STRLENGTH=$(echo -n $test_seq1 | wc -m)

echo ">seq1" > "measuring.fa"
echo ">seq2" > "measuring2.fa"
for i in $(seq 1 $3); do
    echo $test_seq1 >> "measuring.fa"
    echo $test_seq2 >> "measuring2.fa"

    >&2 echo $((STRLENGTH*i)) | tr -d '\n'
    /usr/bin/time -f "%e" ./global_alignment.py measuring.fa measuring2.fa score_matrix -a $1 $2
    echo $time | tr -d '\n'

done
