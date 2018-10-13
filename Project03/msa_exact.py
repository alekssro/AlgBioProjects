#!/usr/bin/env python3
# Usage: ./msa_exact.py score_matrix 5 test.fa
# Requires Alignment.py in the directory

from Alignment import Alignment
from Alignment import GetArguments
import sys

arguments = GetArguments(sys.argv)
test = Alignment(arguments.seqs, arguments.score_matrix, arguments.gapcost)
score = test.sp_exact_3()
alignm = test.backtrack_msa_exact()

for i in range(len(alignm)):
    print(">", arguments.heads[i])
    print(test.num_to_sequence(alignm[i]), "\n")
