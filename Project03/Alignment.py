import sys
from Bio import SeqIO

# We define a class for the data in the exercise with the methods we will need, in order to access easily to them
class Alignment:

    N = []  # substitution matrix (given)
    gap_cost = 100
    listSeq = []
    T = []
    seqs_Nums = []
    alignedSeqs = []
    score = 0

    def __init__(self, listSeq, N, gap_cost, alignType="-l"):

        self.listSeq = listSeq
        self.N = N
        self.gap_cost = gap_cost
        self.gap_params = [0, self.gap_cost]
        self.alignType = alignType

        for seq in self.listSeq:
            self.seqs_Nums.append(self.sequence_to_num(seq))

    def sequence_to_num(self, sequence):
        sequence = sequence.upper()
        sequence = sequence.replace('A', '0')
        sequence = sequence.replace('C', '1')
        sequence = sequence.replace('G', '2')
        sequence = sequence.replace('T', '3')

        return sequence

    def num_to_sequence(self, number):
        number = number.upper()
        number = number.replace('0', 'A')
        number = number.replace('1', 'C')
        number = number.replace('2', 'G')
        number = number.replace('3', 'T')

        return number

    def initMatrix(self, rows, columns, value=float("inf")):
        return [[value for i in range(columns)] for j in range(rows)]

    def initMatrix3d(self, n1, n2, n3, value=float("inf")):
        return [[[value for k in range(n3)] for j in range(n2)] for i in range(n1)]

    def gapcost(self, k):     # k: gap length
        # check type of gap cost depending on self.alignType value; return corresponding gap cost value
        if self.alignType == "-c":
            return self.gap_params[0]

        elif self.alignType == "-a" or self.alignType == "-l":
            return self.gap_params[0] * k + self.gap_params[1]

    def subCost(self, X, Y, i, j):
        return self.N[int(X[i-1])][int(Y[j-1])]

    def calcD(self, i, j):
        d1 = d2 = float("inf")
        d1 = self.T[i-1][j] + self.gapcost(1)
        if (i > 1) and (j >= 0):
            d2 = self.D[i-1][j] + self.gapcost(0)
        self.D[i][j] = min(d1, d2)
        return self.D[i][j]

    def calcI(self, i, j):
        i1 = i2 = float("inf")
        i1 = self.T[i][j-1] + self.gapcost(1)
        if (i >= 0) and (j > 1):
            i2 = self.Ins[i][j-1] + self.gapcost(0)
        self.Ins[i][j] = min(i1, i2)
        return self.Ins[i][j]

    # Returns the sp-score of the MSA stored in the FASTA file 'filename'
    ###########################################################################
    def compute_sp_score(self):

        # Compute the score of each induced pairwise alignment
        score = 0
        for i in range(len(self.M)):
            for j in range(i+1, len(self.M)):
                if len(self.M[i]) != len(self.M[j]):
                    print("ERROR: Rows", i, "and", j, "have different lengths.")
                    sys.exit(1)
                for k in range(len(self.M[i])):
                    print(self.M[i][k], self.M[j][k])
                    score = score + self.N[self.M[i][k]][self.M[j][k]]

        return score

    def findCenterStr(self):
        sumOfScores = float("inf")
        centerStr = 0
        scores = self.initMatrix(len(self.seqs_Nums), len(self.seqs_Nums), 0)
        for i in range(len(self.seqs_Nums)):
            for j in range(len(self.seqs_Nums)):
                if i != j:
                    scores[i][j] = self.affine_align(self.seqs_Nums[i], self.seqs_Nums[j])

            if sum(scores[i]) < sumOfScores:
                sumOfScores = sum(scores[i])
                centerStr = i

        return centerStr

    def extendM(self, optAlign):

        if not self.M:
            self.M = optAlign   # first optimal alignment = M (initialize)
            return

        i = 0
        while self.M[0] != optAlign[0]:

            if i < len(self.M[0]):

                if self.M[0][i] == "-":
                    optAlign[0] = optAlign[0][:i] + '-' + optAlign[0][i:]
                    optAlign[1] = optAlign[1][:i] + '-' + optAlign[1][i:]
                elif optAlign[0][i] == "-":
                    for j in range(len(self.M)):
                        self.M[j] = self.M[j][:i] + "-" + self.M[j][i:]

                i += 1
            elif optAlign[0][i] == "-":
                for j in range(len(self.M)):
                    self.M[j] = self.M[j][:i] + "-" + self.M[j][i:]

        self.M.append(optAlign[1])

    def multiple_align(self):
        self.centerStr = self.findCenterStr()       # index for the center string
        self.M = None
        # self.M = self.initMatrix(2, len(self.seqs_Nums[self.centerStr]))
        optAlign = []
        for i in range(len(self.seqs_Nums)):
            if i == self.centerStr:
                continue

            self.affine_align(self.seqs_Nums[self.centerStr], self.seqs_Nums[i])
            optAlign = self.backtrack_iterative()
            self.extendM(optAlign)

        return(self.M)

    def sp_exact_3(self):

        A = self.seqs_Nums[0]
        B = self.seqs_Nums[1]
        C = self.seqs_Nums[2]
        self.T = self.initMatrix3d(len(A)+1, len(B)+1, len(C)+1)

        for i in range(len(A)+1):
            for j in range(len(B)+1):
                for k in range(len(C)+1):
                    v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float("inf")
                    if i==0 and j==0 and k==0:
                        v0 = 0
                    if i>0 and j>0 and k>0:
                        v1 = self.T[i-1][j-1][k-1] + self.subCost(A, B, i, j) + \
                            self.subCost(B, C, j, k) + self.subCost(A, C, i, k)
                    if i>0 and j>0 and k>=0:
                        v2 = self.T[i-1][j-1][k] + self.subCost(A, B, i, j) + self.gapcost(0) + self.gapcost(0)
                    if i>0 and j>=0 and k>0:
                        v3 = self.T[i-1][j][k-1] + self.gapcost(0) + self.subCost(A, C, i, k) + self.gapcost(0)
                    if i>=0 and j>0 and k>0:
                        v4 = self.T[i][j-1][k-1] + self.gapcost(0) + self.gapcost(0) + self.subCost(B, C, j, k)
                    if i>0 and j>=0 and k>=0:
                        v5 = self.T[i-1][j][k] + self.gapcost(0) + self.gapcost(0)
                    if i>=0 and j>0 and k>=0:
                        v6 = self.T[i][j-1][k] + self.gapcost(0) + self.gapcost(0)
                    if i>=0 and j>=0 and k>0:
                        v7 = self.T[i][j][k-1] + self.gapcost(0) + self.gapcost(0)

                    self.T[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)

        return self.T[i][j][k]

    def backtrack_msa_exact(self):

        A = self.seqs_Nums[0]
        B = self.seqs_Nums[1]
        C = self.seqs_Nums[2]
        a = b = c = ''
        i = len(self.T) - 1
        j = len(self.T[0]) - 1
        k = len(self.T[0][0]) - 1
        while (i > 0 or j > 0 or k > 0):

            if (i > 0 and j > 0 and k > 0) and (self.T[i][j][k] == self.T[i-1][j-1][k-1] +
                    self.subCost(A, B, i, j) + self.subCost(B, C, j, k) + self.subCost(A, C, i, k)):
                # No GAPS
                a = A[i-1] + a
                b = B[j-1] + b
                c = C[k-1] + c
                i -= 1
                j -= 1
                k -= 1
            elif (i > 0 and j > 0 and k >= 0) and (self.T[i][j][k] == self.T[i-1][j-1][k] +
                    self.subCost(A, B, i, j) + self.gapcost(0) + self.gapcost(0)):
                # Gap in C
                a = A[i-1] + a
                b = B[j-1] + b
                c = ("-") + c
                i -= 1
                j -= 1

            elif (i>0 and j>=0 and k>0) and (self.T[i][j][k] == self.T[i-1][j][k-1] + self.gapcost(0) +
                    self.subCost(A, C, i, k) + self.gapcost(0)):
                # Gap in B
                a = A[i-1] + a
                b = ("-") + b
                c = C[k-1] + c
                i -= 1
                k -= 1

            elif (i>=0 and j>0 and k>0) and (self.T[i][j][k] == self.T[i][j-1][k-1] + self.gapcost(0) +
                    self.gapcost(0) + self.subCost(B, C, j, k)):
                # Gap in A
                a = ("-") + a
                b = B[j-1] + b
                c = C[k-1] + c
                j -= 1
                k -= 1

            elif(i > 0 and j >= 0 and k >= 0) and (self.T[i][j][k] == self.T[i-1][j][k] +
                    self.gapcost(0) + self.gapcost(0)):
                # Gap in B and C
                a = A[i-1] + a
                b = ("-") + b
                c = ("-") + c
                i -= 1

            elif(i>=0 and j>0 and k>=0) and (self.T[i][j][k] == self.T[i][j-1][k] +
                    self.gapcost(0) + self.gapcost(0)):

                # Gap in A and C
                a = ("-") + a
                b = B[j-1] + b
                c = ("-") + c
                j -= 1

            elif(i >= 0 and j >= 0 and k > 0) and (self.T[i][j][k] == self.T[i][j][k-1] +
                    self.gapcost(0) + self.gapcost(0)):
                # Gap in A and B
                a = ("-") + a
                b = ("-") + b
                c = C[k-1] + c
                k -= 1

        return a, b, c

    def affine_align(self, seq1, seq2):
        """Global alignment with affine or linear penalties. We assume we are minimizing."""
        self.D = self.initMatrix(len(seq1) + 1, len(seq2) + 1)
        self.Ins = self.initMatrix(len(seq1) + 1, len(seq2) + 1)
        self.T = self.initMatrix(rows=len(seq1) + 1, columns=len(seq2) + 1)
        self.A, self.B = seq1, seq2

        self.T[0][0] = 0
        # Rest of the matrix
        for j in range(0, len(self.T[0])):
            for i in range(1, len(self.T)):
                v1 = v2 = v3 = float("inf")

                if (i > 0) and (j > 0):     # Diagonal arrow
                    v1 = self.T[i-1][j-1] + self.N[int(self.A[i-1])][int(self.B[j-1])]
                if (i > 0) and (j >= 0):    # Vertical arrow
                    v2 = self.calcD(i, j)
                if (i >= 0) and (j > 0):    # Horizontal arrow
                    v3 = self.calcI(i, j)

                self.T[i][j] = min(v1, v2, v3)

        self.score = self.T[i][j]

        return self.T[i][j]

    def backtrack_iterative(self):

        a = b = ''
        i = len(self.T) - 1
        j = len(self.T[0]) - 1
        while (i>0 or j>0):
            if (i > 0 and j > 0) and (self.T[i][j] == self.T[i-1][j-1] + self.N[int(self.A[i-1])][int(self.B[j-1])]):
                # optimal alignment of A[1..i] and B[1..j] ends in a sub-column
                a = self.A[i-1] + a
                b = self.B[j-1] + b
                i -= 1
                j -= 1
            else:
                # optimal alignment of A[1..i] and B[1..j] ends in a del- or in-block
                k = 1
                while True:
                    if (i >= k) and self.T[i][j] == self.T[i-k][j] + self.gapcost(k):
                        # optimal alignment of A[1..i] and B[1..j] ends in del-block of length k
                        # “output columns (A[i] A[i-1] … A[i-k+1], -- … --)”
                        a = self.A[i-k:i] + a
                        b = ("-"*k) + b
                        i = i - k
                        break
                    elif (j >= k) and self.T[i][j] == self.T[i][j-k] + self.gapcost(k):
                        # optimal alignment of A[1..i] and B[1..j] ends in a in-block of length k
                        # “output columns (-- … --, B[j]B[j-1] … B[j-k+1])”
                        a = ("-"*k) + a
                        b = self.B[j-k:j] + b
                        j = j - k
                        break
                    else:
                        k = k + 1

        return [a, b]

    def align(self):

        self.score = self.sp_exact_3()
        self.a, self.b, self.c = self.backtrack_msa_exact()
        # self.backtrack_iterative()
        self.a = self.num_to_sequence(self.a)
        self.b = self.num_to_sequence(self.b)
        self.c = self.num_to_sequence(self.c)

class GetArguments:

    argList = []
    seqs = []  # make it generic for next projects ??
    score_matrix = []
    gap_params = 100
    n_chr = 0       # number of different characters in the score matrix
    characters = []

    def __init__(self, argList):

        self.output = False
        if argList[len(argList)-1] == "-o":
            self.output = True
        self.argList = argList  # get class input (list of arguments provided)

        self.score_matrix = self.readScoreMatrix(argList[1])
        self.gapcost = float(argList[2])
        self.seqs = self.read_fasta(argList[3])

    def check_num_args(self):
        # Check correct number of arguments
        if (len(self.argList) < 5):
            print("Error: Insufficient number of arguments.")
            # print("\tGiven", len(self.argList)-1, "when expecting at least 4.")
            print("\tType '", self.argList[0], "-h' for help.")
            sys.exit(1)
        else:
            # print "Number of arguments correct (>4)."
            pass

    def check_alignment_type(self, alignType):
        # Check if alignment type is allowed and return corresponding flag to be used in Alignment

        if len(self.gap_params) > 0:
            if (alignType == "-c" or alignType == "-constant"):
                # print("Performing alignment with constant gap cost ->", self.gap_params[0])
                return "-c"
            elif (alignType == "-l" or alignType == "-linear") and len(self.gap_params) > 0:
                self.gap_params = [0, self.gap_params[0]]
                # print("Performing alignment with linear gap cost -> b =", self.gap_params[1])
                return "-l"
            elif (alignType == "-a" or alignType == "-affine"):
                if len(self.gap_params) > 1:
                    # print("Performing alignment with affine gap cost -> extend gap cost =", self.gap_params[0],
                            # "\tstart gap block =", self.gap_params[1])
                    return "-a"
        else:
            print("Error: Incorrect type of alignment or insufficient number of parameters to perform the alignment.")
            print("\tType '", self.argList[0], "-h' for help.")
            sys.exit(1)

    def read_fasta(self, infile):
        seqs = []
        self.heads = []
        for record in SeqIO.parse(infile, "fasta"):
            self.heads.append(str(record.id))
            seqs.append(str(record.seq))

        return seqs

    def readScoreMatrix(self, infile):

        lines = [line.rstrip('\n') for line in open(infile)]    # lines to elements in list
        words = []
        for line in lines:
            words.append(line.split())          # devide the lines into characters

        self.n_chr = int(words[0][0])

        # initialize score matrix to fill
        score_matrix = [[0 for i in range(self.n_chr)] for j in range(self.n_chr)]  # inizilize score matrix

        # safe different characters in self.characters and score_matrix
        for i in range(1, self.n_chr+1):
            self.characters.append(words[i][0])
            for j in range(1, self.n_chr+1):
                score_matrix[i-1][j-1] = int(words[i][j])

        return score_matrix
