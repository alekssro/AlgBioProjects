import sys
from Bio import SeqIO

# We define a class for the data in the exercise with the methods we will need, in order to access easily to them
class Alignment:

    N = []  # substitution matrix (given)
    alignType = ""
    gap_params = []
    k = 0       # gap length (0 for constant so default)
    listSeq = []
    seq1 = ""
    seq2 = ""
    T = []
    A, B = "", ""  # A and B are sequences translated to numbers (0-3)
    a, b = "", ""  # a and b are the sequences after the alignment
    score = 0

    def __init__(self, listSeq, N, alignType, gap_params):

        self.listSeq = listSeq
        self.seq1 = listSeq[0]
        self.seq2 = listSeq[1]
        self.N = N
        self.alignType = alignType
        self.k = 0      # gap length (0 for constant so default; can change inside cost function)
        self.gap_params = gap_params
        self.T = self.initMatrix(rows=len(self.seq1) + 1, columns=len(self.seq2) + 1, value=float("inf"))  # inizilize score matrix

        self.A, self.B = self.sequence_to_num(self.seq1), self.sequence_to_num(self.seq2)

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

    def initMatrix(self, rows, columns, value=float("-inf")):
        return [[value for i in range(columns)] for j in range(rows)]

    def gapcost(self, k):     # k: gap length
        # check type of gap cost depending on self.alignType value; return corresponding gap cost value
        if self.alignType == "-c":
            return self.gap_params[0]

        elif self.alignType == "-a" or self.alignType == "-l":
            return self.gap_params[0] * k + self.gap_params[1]

    def affine_align(self):
        """Global alignment with affine or linear penalties. We assume we are minimizing."""
        self.D = self.initMatrix(len(self.seq1) + 1, len(self.seq2) + 1, value=0)
        self.Ins = self.initMatrix(len(self.seq1) + 1, len(self.seq2) + 1, value=0)

        # Rest of the matrix
        for j in range(0, len(self.T[0])):
            v1 = v2 = v3 = float("inf")
            for i in range(0, len(self.T)):

                if (i == 0) and (j == 0):
                    self.T[i][j] = 0
                    break
                if (i > 0) and (j > 0):     # Diagonal arrow
                    v1 = self.T[i-1][j-1] + self.N[int(self.A[i-1])][int(self.B[j-1])]
                if (i > 0) and (j >= 0):    # Vertical arrow
                    d1 = self.T[i-1][j] + self.gapcost(1)
                    d2 = float("inf")
                    if (i > 1) and (j >= 0):
                        d2 = self.D[i-1][j] + self.gapcost(0)
                    self.D[i][j] = min(d1, d2)
                    v2 = self.D[i][j]
                if (i >= 0) and (j > 0):    # Horizontal arrow
                    i1 = self.T[i][j-1] + self.gapcost(1)
                    i2 = float("inf")
                    if (i >= 0) and (j > 1):
                        i2 = self.Ins[i][j-1] + self.gapcost(0)
                    self.Ins[i][j] = min(i1, i2)
                    v3 = self.Ins[i][j]

                self.T[i][j] = min(v1, v2, v3)

        self.score = self.T[i][j]

        return self.T[i][j]

    def RecurBackTrack(self, i, j, k=-5):

        if(i > 0) and (j > 0) and self.T[i][j] == (self.T[i-1][j-1] + self.N[int(self.A[i-1])][int(self.B[j-1])]):
            self.a = self.A[i-1] + self.a
            self.b = self.B[j-1] + self.b
            self.RecurBackTrack(i-1, j-1)
        elif(i > 0) and (j >= 0) and (self.T[i][j] == self.T[i-1][j] + self.gapcost(1) or self.T[i][j] == self.D[i-1][j] + self.gapcost(0)):
            self.a = self.A[i-1] + self.a
            self.b = "-" + self.b
            self.RecurBackTrack(i-1, j)
        elif(i >= 0) and (j > 0) and (self.T[i][j] == self.T[i][j-1] + self.gapcost(1) or self.T[i][j] == self.Ins[i][j-1] + self.gapcost(0)):
            self.a = "-" + self.a
            self.b = self.B[j-1] + self.b
            self.RecurBackTrack(i, j-1)

        return self.a, self.b

    def align(self):

        self.affine_align()

        # self.iterBackTrack()
        # self.a = self.num_to_sequence(self.a)
        # self.b = self.num_to_sequence(self.b)
        # print(self.a, "\n", self.b)
        # self.a, self.b = "", ""
        self.RecurBackTrack(len(self.seq1), len(self.seq2))
        self.a = self.num_to_sequence(self.a)
        self.b = self.num_to_sequence(self.b)

class GetArguments:

    argList = []
    seq1, seq2 = '', ''  # make it generic for next projects ??
    score_matrix = []
    alignment_type = ""
    gap_params = []
    n_chr = 0       # number of different characters in the score matrix
    characters = []

    def __init__(self, argList):

        self.output = False
        if argList[len(argList)-1] == "-o":
            self.output = True
        self.argList = argList  # get class input (list of arguments provided)
        self.check_num_args()        # will exit program if arguments provided are insufficient
        self.seq1 = self.read_fasta(argList[1])
        self.seq2 = self.read_fasta(argList[2])
        self.score_matrix = self.readScoreMatrix(argList[3])

        if self.output:
            end_gap_params = len(argList) - 1
        else:
            end_gap_params = len(argList)

        self.gap_params = [float(x) for x in argList[5:end_gap_params]]   # Get gap cost arguments as float
        self.alignment_type = self.check_alignment_type(argList[4])

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
        fasta_record = SeqIO.read(infile, 'fasta')
        sequence = str(fasta_record.seq)
        return sequence

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
