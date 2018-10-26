#!/usr/bin/env python3
####################################################################################################
# Program used to calculate the RF distance between two unrooted evolutionary trees
# over the same set of species. We use Day's algorithm to calculate the distance.
#
# Input: takes two evolutionary trees in Newick format (also referred to as 'New Hampshire format')
# Output: RF distance between the two inputed trees
#
# Usage: ./rfdist.py tree1.new tree2.new
#
# Arguments:
#   - tree1.new, tree2.new: files containing unrooted trees in Newick format.
####################################################################################################

# Needed libraries
from FiloTree import GetArguments
from FiloTree import FiloTree
import sys

# Read arguments:
arguments = GetArguments(sys.argv)

tree1 = arguments.trees1
tree2 = arguments.trees2

# Get distance
Trees = FiloTree(tree1, tree2)

Trees.numerDF()

Trees.dist()
