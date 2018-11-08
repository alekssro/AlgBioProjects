#!/bin/bash
# Script to measure the running time of our implementation of global_alignment.py

last=10
for (( i = 0; i < 38; i++ )); do

    head -n $last Data/Alignments/Normal/clustalo_patbase_aitbas.stockholm > Timing/test_alignm.stockholm
    echo "" >> Timing/test_alignm.stockholm
    echo '//' >> Timing/test_alignm.stockholm
    (./Programs/quicktree_1.1/bin/quicktree Timing/test_alignm.stockholm > "Timing/Trees/test_tree_quicktree.newick") 2>trash
    (./Programs/rapidNJ/bin/rapidnj -i sth Timing/test_alignm.stockholm | sed -e "s/'//g" > "Timing/Trees/test_tree_rapidNJ.newick") 2> trash

    >&2 echo $((last-2)) | tr -d '\n'
    (/usr/bin/time -f "%e" ./rfdist.py Timing/Trees/test_tree_quicktree.newick Timing/Trees/test_tree_rapidNJ.newick)

    last=$((last+10))
done
