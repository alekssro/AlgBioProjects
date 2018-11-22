#!/bin/bash
# Script to experiment

for dist_Tree in ../Data/distance_matrices/unique_distance_matrices/*
do
    filename=$(basename -- "$dist_Tree")
    # echo $filename

    # time quicktree
    time_QT="$( TIMEFORMAT='%3R';time (../Software/quicktree_1.1/bin/quicktree $dist_Tree > Trees/${filename}_quicktree.newick) 2>&1 1>/dev/null)"

    # time Neighbor-joining
    start=$(($(date +%s%N)/1000000))
    (../Software/rapidNJ/bin/rapidnj -i pd $dist_Tree -v | sed -e "s/'//g" > Trees/${filename}_rapinj.newick) 2>/dev/null
    end=$(($(date +%s%N)/1000000))
    time_RNJ=$(((end-start)))

    # time our Program
    time_emar="$( TIMEFORMAT='%3R';time (../emar-nj.py $dist_Tree > Trees/${filename}_emar.newick) 2>&1 1>/dev/null)"
    echo "${filename%%_*}" $time_QT $time_RNJ $ time_emar
done
