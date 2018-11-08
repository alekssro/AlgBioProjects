#!/bin/bash
# Script to create trees using quicktree and rapidNJ from the multiple alignments in Data/Alignments/

for tree1 in Data/Trees/Permuted/*
do
    array=()
    for tree2 in Data/Trees/Permuted/*
    do
        echo $(./rfdist.py $tree1 $tree2) | tr -d '\n'
        echo " " | tr -d '\n'
    done
    echo ""
done
