#!/bin/bash
# Script to evaluate our implementation of global_alignment.py to answer questions 3 and 4 from project2_eval.txt


for seq1 in sequences/*
do
    array=()
    for seq2 in sequences/*
    do
        echo $(./global_alignment.py $seq1 $seq2 score_matrix -a $1 $2) | tr -d '\n'
        echo " " | tr -d '\n'
    done
    echo ""
done
