# from Bio import Phylo


class NJtree(object):
    """Methods for constructing a NJ tree"""

    dist_matrix = []
    taxa_names = []
    N = []
    tree = ''
    clades = {}

    def __init__(self, dist_matrix, taxa_names):

        self.dist_matrix = dist_matrix
        self.taxa_names = taxa_names

        self.buildNJ()

    def buildNJ(self):

        while len(self.dist_matrix) > 3:

            # Get closest pair
            pair = self.getClosestPair()
            i, j = pair

            self.addNode(i, j)
            self.updateDist_matrix(i, j)

            name_j = self.taxa_names.pop(j)
            name_i = self.taxa_names.pop(i)

            self.taxa_names.append(name_i + name_j)

        self.addFinalNode()

    def addFinalNode(self):

        v_i = 1 / 2 * (self.dist_matrix[0][1] + self.dist_matrix[0][2] - self.dist_matrix[1][2])
        v_j = 1 / 2 * (self.dist_matrix[1][0] + self.dist_matrix[1][2] - self.dist_matrix[0][2])
        v_m = 1 / 2 * (self.dist_matrix[2][0] + self.dist_matrix[2][1] - self.dist_matrix[0][1])

        if self.taxa_names[0] in self.clades:
            clade_1 = self.clades[self.taxa_names[0]] + ":" + str(v_i)
        else:
            clade_1 = self.taxa_names[0] + ":" + str(v_i)

        if self.taxa_names[1] in self.clades:
            clade_2 = self.clades[self.taxa_names[1]] + ":" + str(v_j)
        else:
            clade_2 = self.taxa_names[1] + ":" + str(v_j)

        if self.taxa_names[2] in self.clades:
            clade_3 = self.clades[self.taxa_names[2]] + ":" + str(v_m)
        else:
            clade_3 = self.taxa_names[2] + ":" + str(v_m)

        self.tree = "(" + clade_1 + "," + clade_2 + "," + clade_3 + ")" + ";"

    def addNode(self, i, j):

        w_i = 1 / 2 * (self.dist_matrix[i][j] + self.calc_ri(i) - self.calc_ri(j))
        w_j = 1 / 2 * (self.dist_matrix[i][j] + self.calc_ri(j) - self.calc_ri(i))

        leave_1 = self.taxa_names[i] + ":" + str(w_i)
        leave_2 = self.taxa_names[j] + ":" + str(w_j)

        if self.taxa_names[i] in self.clades and self.taxa_names[j] in self.clades:

            self.clades[self.taxa_names[i] + self.taxa_names[j]] = "(" + self.clades[self.taxa_names[i]] + ":" + str(
                w_i) + "," + self.clades[self.taxa_names[j]] + "):" + str(w_j)
            del self.clades[self.taxa_names[i]]
            del self.clades[self.taxa_names[j]]

        elif self.taxa_names[i] in self.clades:

            self.clades[self.taxa_names[i] + self.taxa_names[j]] = "(" + self.clades[self.taxa_names[i]] + ":" + str(
                w_i) + "," + leave_2 + ")"
            del self.clades[self.taxa_names[i]]

        elif self.taxa_names[j] in self.clades:

            self.clades[self.taxa_names[i] + self.taxa_names[j]] = "(" + self.clades[self.taxa_names[j]] + ":" + str(
                w_j) + "," + leave_1 + ")"
            del self.clades[self.taxa_names[j]]

        else:
            self.clades[self.taxa_names[i] + self.taxa_names[j]] = "(" + leave_1 + "," + leave_2 + ")"

    def updateDist_matrix(self, i, j):

        # Calculate new row/column distance
        d_k = []
        for m in range(len(self.dist_matrix)):
            if m != i and m != j:
                dist = 1 / 2 * \
                    (self.dist_matrix[i][m] + self.dist_matrix[j]
                     [m] - self.dist_matrix[i][j])
                d_k.append(dist)

        # Remove i and j rows/columns
        # Remove rows
        del self.dist_matrix[max(i, j)]
        del self.dist_matrix[min(i, j)]

        # Remove columns
        for m in range(len(self.dist_matrix)):
            del self.dist_matrix[m][max(i, j)]
            del self.dist_matrix[m][min(i, j)]
            self.dist_matrix[m].append(d_k[m])

        # Add new row
        d_k.append(0)
        self.dist_matrix.append(d_k)

    def getClosestPair(self):

        # inizilize N matrix
        N = [[0 for column in range(len(self.dist_matrix))]
             for row in range(len(self.dist_matrix[0]))]
        min_val = float("inf")
        min_indexes = [0, 0]

        for i in range(len(N)):
            for j in range(len(N[i])):
                N[i][j] = self.dist_matrix[i][j] - \
                    (self.calc_ri(i) + self.calc_ri(j))

                if i != j and N[i][j] < min_val:
                    min_val = N[i][j]
                    min_indexes = [i, j]

        return min_indexes

    def calc_ri(self, i):

        r = (1 / (len(self.dist_matrix) - 2)) * sum(self.dist_matrix[i])

        return r


class ReadPhylip:
    """Extracts information from a phylip format file
    Input: Phylip-format file name
    Places the matrix values in 'phylip_matrix' and the names in 'characters'"""

    argList = []
    phylip_matrix = []
    n_chr = 0       # number of different characters in the score matrix
    characters = []

    def __init__(self, argList):

        # self.output = False
        # if argList[len(argList) - 1] == "-o":
        #     self.output = True
        self.argList = argList  # get class input (list of arguments provided)

        self.phylip_matrix = self.readPhylipMatrix(argList[1])

    def readPhylipMatrix(self, infile):

        lines = [line.rstrip('\n') for line in open(
            infile)]    # lines to elements in list
        words = []
        for line in lines:
            # devide the lines into characters
            words.append(line.split())

        self.n_chr = int(words[0][0])

        # initialize score matrix to fill
        phylip_matrix = [[0 for i in range(self.n_chr)] for j in range(
            self.n_chr)]  # inizilize score matrix

        # safe different characters in self.characters and phylip_matrix
        for i in range(1, self.n_chr + 1):
            self.characters.append(words[i][0])
            for j in range(1, self.n_chr + 1):
                phylip_matrix[i - 1][j - 1] = float(words[i][j])

        return phylip_matrix
