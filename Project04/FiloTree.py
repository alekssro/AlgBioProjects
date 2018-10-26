from Bio import Phylo
import sys


class FiloTree:

    tree1 = None
    tree2 = None
    distance = float("Inf")
    leaf2num = {}

    def __init__(self, tree1, tree2):
        self.tree1 = tree1
        self.tree2 = tree2

        # Phylo.draw_ascii(self.tree1)
        # Phylo.draw_ascii(self.tree2)

        self.distance = self.dist()

    def dist(self):
        self.rootTrees()
        self.numerDF()
        self.getIntervals()

    def rootTrees(self):
        for clade in self.tree1.find_clades():
            if clade.name:
                root = clade.name   # save first leaf as root
                break
        self.tree1.root_with_outgroup(root)     # root tree1
        self.tree2.root_with_outgroup(root)     # root tree2
        print(root)

    def numerDF(self):
        # method for numbering each leaf in T1; Depth-First search
        num = 0
        for clade in self.tree1.find_clades():
            if clade.name:
                num += 1
                self.leaf2num[clade.name] = num

    # TODO: method for annotating internal nodes with interval (numbering)
    #       if “max – min + 1 = size” then it is an interval, potential candidate
    def getIntervals(self):
        pass

    # TODO: method to identify shared intervals
    def sort_count(self, arg):
        pass


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
