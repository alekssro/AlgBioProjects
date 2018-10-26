from Bio import Phylo
import sys

class FiloTree:

    tree1 = None
    tree2 = None
    distance = float("Inf")

    def __init__(self, tree1, tree2):
        self.tree1 = tree1
        self.tree2 = tree2

        # Phylo.draw_ascii(self.tree1)
        # Phylo.draw_ascii(self.tree2)

        self.distance = self.dist()

    def dist(self):
        pass

    # TODO: method for rooting a tree at a defined leaf
    # TODO: method for numbering each leaf in T1 Depth-First
    # TODO: method for associating leaves in T2 with the numbering done in T1
    # TODO: method for annotating internal nodes with interval (numbering)
    #       if “max – min + 1 = size” then it is an interval, potential candidate
    # TODO: method to identify shared intervals

class GetArguments:

    def __init__(self, argList):
        self.argList = argList
        self.check_num_args()
        self.trees1 = self.readTree(argList[1])
        self.trees2 = self.readTree(argList[2])

    def check_num_args(self):
        # Check correct number of arguments
        if (len(self.argList) != 3):
            print("Error: Incorrect number of arguments. Expected 2.")
            print("\tUsage:", self.argList[0], "tree1.new tree2.new")
            sys.exit(1)
        else:
            pass

    def readTree(self, treeFile):

        tree = Phylo.read(treeFile, 'newick')

        return tree
