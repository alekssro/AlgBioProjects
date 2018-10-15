from Bio import Phylo


class FiloTree:

    tree1 = None
    tree2 = None
    distance = float("Inf")

    # TODO: Define methods needed for calculating the distance
    def __init__(self, tree1, tree2):
        self.tree1 = tree1
        self.tree2 = tree2

        Phylo.draw_ascii(self.tree1)
        Phylo.draw_ascii(self.tree2)

        self.distance = self.dist()

    def dist(self):
        pass


class GetArguments:

    def __init__(self, argList):
        self.argList = argList
        self.trees1 = self.readTree(argList[1])
        self.trees2 = self.readTree(argList[2])

    def checkArgs(self, arg):
        # TODO: implement checkArgs
        pass

    def readTree(self, treeFile):

        trees = []
        for record in Phylo.parse(treeFile, "newick"):
            trees.append(record)

        return trees
