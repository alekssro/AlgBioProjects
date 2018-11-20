#!/usr/bin/env python3
####################################################################################################
# Program used to create a Neighbor-Joining (NJ) tree from a given distance matrix
#
# Input: takes a distance matrix in Phylip format
# Output: NJ tree in newick format
#
# Usage: ./emar-nj.py distance-matrix.phy
#
# Arguments:
#   - distance-matrix.phy: distance matrix in Phylip format
####################################################################################################

# Needed libraries
from NJ import ReadPhylip
from NJ import NJtree
import sys

# Read arguments:
arguments = ReadPhylip(sys.argv)

dist_matrix = arguments.phylip_matrix
taxa = arguments.characters
print(taxa)

# Build NJ tree
Tree = NJtree(dist_matrix, taxa)

print(Tree)
