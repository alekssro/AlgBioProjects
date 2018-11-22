#!/bin/bash
# Script to experiment
# Execute from Test directory

# Header
echo "N time_QT time_RNJ time_emar QT_emar RNJ_emar QT_RNJ"
for dist_Tree in ../Data/distance_matrices/unique_distance_matrices/*
do
    filename=$(basename -- "$dist_Tree")
    # echo $filename

    # time quicktree
    time_QT="$( TIMEFORMAT='%3R';time (../Software/quicktree_1.1/bin/quicktree $dist_Tree > Trees/quicktree.newick) 2>&1 1>/dev/null)"

    # time Neighbor-joining # Note: it's in miliseconds
    start=$(($(date +%s%N)/1000000))
    (../Software/rapidNJ/bin/rapidnj -i pd $dist_Tree -v | sed -e "s/'//g" > Trees/rapidnj.newick) 2>/dev/null
    end=$(($(date +%s%N)/1000000))
    time_RNJ=$(((end-start)))

    # time our Program
    time_emar="$( TIMEFORMAT='%3R';time (../emar-nj.py $dist_Tree > Trees/emar.newick) 2>&1 1>/dev/null)"

    # Distance between Trees
    QT_emar=$(../Software/rfdist.py Trees/quicktree.newick Trees/emar.newick)
    RNJ_emar=$(../Software/rfdist.py Trees/rapidnj.newick Trees/emar.newick)
    QT_RNJ=$(../Software/rfdist.py Trees/quicktree.newick Trees/rapidnj.newick)

    echo "${filename%%_*}" $time_QT $time_RNJ $time_emar $QT_emar $RNJ_emar $QT_RNJ

done
