#!/usr/bin/env python3
# This project is about implementing and experimenting with pairwise sequence comparison methods
# to compute optimal global alignments of two sequences where the object is to minimize a cost.
# Usage: global_alignment.py seq1.fasta seq2.fasta score_matrix -alignment_type b [a]
# Arguments:
#   - seq1.fasta: fasta file containing sequence 1.
#   - seq2.fasta: fasta file containing sequence 2.
#   - score_matrix: file containing the score matrix used for the alignment. In "Phylip-like" format
#   - -alignment_type: type of alignment to be performed
#       for constant gap cost use -c or -constant
#       for linear gap cost use -l or -linear
#       for affine gap cost use -a or -affine
#   - b, a: parameters for gap cost function
#       b -> constant gap cost or slope when performing linear/affine gap constant (extension penalty)
#       a -> instersect for affine gap cost (opening gap penalty)
#   - -o: output alignment. if missing then outputs optimal score

from Alignment import Alignment
from Alignment import GetArguments
import sys

if (len(sys.argv) > 1) and sys.argv[1] == "-h":
    print("This project is about implementing and experimenting with pairwise sequence comparison methods to compute optimal global alignments of two sequences where the object is to minimize a cost.\n\
    Usage: global_alignment.py seq1.fasta seq2.fasta score_matrix -alignment_type b [a] [-o]\n\
    Arguments:\n\
      - seq1.fasta: fasta file containing sequence 1.\n\
      - seq2.fasta: fasta file containing sequence 2.\n\
      - score_matrix: file containing the score matrix used for the alignment. In 'Phylip-like' format\n\
      - -alignment_type: type of alignment to be performed\n\
          for linear gap cost use -l or -linear\n\
          for affine gap cost use -a or -affine\n\
      - b, a: parameters for gap cost function\n\
          b -> constant gap cost or slope when performing linear/affine gap constant (extension penalty)\n\
          a -> instersect for affine gap cost (opening gap penalty)\n\
      - -o: output alignment. if missing then outputs optimal score")

arguments = GetArguments(sys.argv)      # parsing arguments
# print(arguments.seq2)
# print(arguments.score_matrix)

sequences = [arguments.seq1, arguments.seq2]

substitution_matrix = arguments.score_matrix
gap_params = arguments.gap_params
alignmentType = arguments.alignment_type

my_alignment = Alignment(sequences, substitution_matrix, alignmentType, gap_params)
my_alignment.align()
if arguments.output:
    print(">seq1")
    print(my_alignment.a)
    print(">seq2")
    print(my_alignment.b)
else:
    print(my_alignment.score)

# # Score matrix print
# for row in range(len(T)):
#     print(T[row])
