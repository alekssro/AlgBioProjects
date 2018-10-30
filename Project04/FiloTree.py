from Bio import Phylo
import sys


class FiloTree:

    tree1 = None
    tree2 = None
    distance = float("Inf")
    leaf2num = {}
    intervals1 = []
    intervals2 = []
    non_intervals = 0

    def __init__(self, tree1, tree2):
        self.tree1 = tree1
        self.tree2 = tree2

        # Phylo.draw_ascii(self.tree1)
        # Phylo.draw_ascii(self.tree2)

        self.distance = self.dist()

    def dist(self):
        self.rootTrees()
        self.numerDF()
        self.intervals1 = self.getIntervals(self.tree1)
        self.intervals2 = self.getIntervals(self.tree2)
        return self.sort_and_count()

    def rootTrees(self):
        for clade in self.tree1.find_clades():
            if clade.name:
                root = clade.name   # save first leaf as root
                break
        self.tree1.root_with_outgroup(root)     # root tree1
        self.tree2.root_with_outgroup(root)     # root tree2
        # print(root)

    def numerDF(self):
        # method for numbering each leaf in T1; Depth-First search
        num = 0
        for clade in self.tree1.find_clades():
            if clade.name:
                num += 1
                self.leaf2num[clade.name] = num

    # TODO: method for annotating internal nodes with interval (numbering)
    #       if “max – min + 1 = size” then it is an interval, potential candidate
    def getIntervals(self, tree):
        intervals = []
        for clade in tree.find_clades():
            if not clade.is_terminal():
                terminals = clade.get_terminals()
                terminalsNums = [self.leaf2num[terminal.name] for terminal in terminals]
                mini = min(terminalsNums)
                maxi = max(terminalsNums)
                size = len(terminalsNums)
                if (maxi - mini + 1) == size:
                    intervals.append([mini, maxi])
                else:
                    self.non_intervals += 1

        return intervals

    def sort_and_count(self):
        shared = 0
        all = self.intervals1 + self.intervals2
        all_sorted = sorted(all)
        for i in range(len(all_sorted)-1):
            if all_sorted[i] == all_sorted[i+1]:
                shared += 1

        print(len(all_sorted), shared, self.non_intervals)
        d = len(all_sorted) - 2 * shared + self.non_intervals
        return d


class GetArguments:

    def __init__(self, argList):
        self.argList = argList
        self.check_num_args()
        self.tree1 = self.readTree(argList[1])
        self.tree2 = self.readTree(argList[2])

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
