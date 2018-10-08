#!/usr/bin/env python3
# Usage: ./msa_approx.py score_matrix 5 test.fa

from Alignment import Alignment
from Alignment import GetArguments
import sys

arguments = GetArguments(sys.argv)
test = Alignment(arguments.seqs, arguments.score_matrix, arguments.gapcost)
alignm = test.multiple_align()

for i in range(len(alignm)):
    print("> ", arguments.heads[test.seqOrder[i]], "\n", test.num_to_sequence(alignm[i]), "\n", sep="")
