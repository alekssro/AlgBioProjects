#!/bin/bash
# Script to create trees using quicktree and rapidNJ from the multiple alignments in Data/Alignments/

trees_normal=(Data/Trees/Normal/*)
trees_permuted=(Data/Trees/Permuted/*)
for (( i=0; i<8; i++ ));
do
    echo $(./rfdist.py ${trees_normal[i]} ${trees_permuted[i]}) | tr -d '\n'
    echo " " | tr -d '\n'
done

echo ""
