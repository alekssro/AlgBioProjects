#!/usr/bin/env python3

from Alignment import Alignment
from Alignment import GetArguments
import sys

arguments = GetArguments(sys.argv)
test = Alignment(arguments.seqs, arguments.score_matrix, arguments.gapcost)
print(test.sp_exact_3())
print(test.multiple_align())
# for x in range(len(test.T)):
#     print(test.T[x])
# print(test.T)
# for l in test.T:
#     print([e[0] for e in l])
