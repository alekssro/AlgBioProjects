#!/bin/bash
# Script to create trees using quicktree and rapidNJ from the multiple alignments in Data/Alignments/

for alignm in Data/Alignments/*
do
    filename=$(basename -- "$alignm")
    filename="${filename%.*}"

    ./Programs/quicktree_1.1/bin/quicktree $alignm > "Data/Trees/${filename}_quicktree.newick"
done

for alignm in Data/Alignments/*
do
    filename=$(basename -- "$alignm")
    filename="${filename%.*}"

    ./Programs/rapidNJ/bin/rapidnj -i sth $alignm | sed -e "s/'//g" > "Data/Trees/${filename}_rapidNJ.newick"
done
