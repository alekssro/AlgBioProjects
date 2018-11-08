#!/bin/bash
# Script to calculate rf-distance between all the combinations in Data/Trees/

for tree1 in Data/Trees/Normal/*
do
    array=()
    for tree2 in Data/Trees/Normal/*
    do
        echo $(./rfdist.py $tree1 $tree2) | tr -d '\n'
        echo " " | tr -d '\n'
    done
    echo ""
done
