#!/usr/bin/env python3
# Usage: ./sp_exact_3.py score_matrix 5 test.fa


from Alignment import Alignment
from Alignment import GetArguments
import sys

arguments = GetArguments(sys.argv)
test = Alignment(arguments.seqs, arguments.score_matrix, arguments.gapcost)
score = test.sp_exact_3()
alignm = test.backtrack_msa_exact()

for i in range(len(alignm)):
    print("> ", arguments.heads[i], "\n", test.num_to_sequence(alignm[i]), "\n", sep="")

# for x in range(len(test.T)):
#     print(test.T[x])
# print(test.T)
# for l in test.T:
#     print([e[0] for e in l])
