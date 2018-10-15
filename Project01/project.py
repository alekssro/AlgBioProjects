# We define a class for the data in the exercise with the methods we will need, in order to access easily to them
class Alignment:

    N = [[10,2,5,2],
        [2,10,2,5],
        [5,2,10,2],
        [2,5,2,10]] # substitution maxtrix (given)
    k = -5      # gap cost (given)
    listSeq = []
    seq1 = ""
    seq2 = ""
    T = []
    A,B = "","" # A and B are sequences translated to numbers (0-3)
    a,b = "","" # a and b are the sequences after the alignment
    score = -1

    def __init__(self, listSeq):

        self.listSeq = listSeq
        self.seq1 = listSeq[0]
        self.seq2 = listSeq[1]
        self.T = [ [ float("-inf") for i in range(len(self.seq2)+1) ] for j in range(len(self.seq1)+1) ] #inizilize score matrix

        self.A,self.B = self.sequence_to_num(self.seq1), self.sequence_to_num(self.seq2)

    def sequence_to_num(self,sequence):
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

    def cost(self,i,j,k=-5):

        if self.T[i][j] == float("-inf"):
            v1 = v2 = v3 = v4 = float("-inf")

            if (i > 0) and (j > 0):
                v1 = self.cost(i-1, j-1, k) + self.N[int(self.A[i-1])][int(self.B[j-1])]
            if (i > 0) and (j >= 0):
                v2 = self.cost(i-1, j, k) + k
            if (i >= 0) and (j > 0):
                v3 = self.cost(i, j-1, k) + k
            if (i == 0) and (j == 0):
                v4 = 0

            self.T[i][j] = max(v1,v2,v3,v4)

        self.score = self.T[i][j]
        return self.T[i][j]

    def RecurBackTrack(self,i, j, k=-5):

        if(i > 0) and (j > 0) and self.T[i][j] == (self.T[i-1][j-1] + self.N[int(self.A[i-1])][int(self.B[j-1])]):
            self.a = self.A[i-1] + self.a
            self.b = self.B[j-1] + self.b
            self.RecurBackTrack(i-1, j-1)
        elif(i > 0) and (j >= 0) and self.T[i][j] == self.T[i-1][j] + k:
            self.a = self.A[i-1] + self.a
            self.b = "-" + self.b
            self.RecurBackTrack(i-1, j)
        elif(i >= 0) and (j > 0) and self.T[i][j] == self.T[i][j-1] + k:
            self.a = "-" + self.a
            self.b = self.B[j-1] + self.b
            self.RecurBackTrack(i, j-1)

        # return self.a, self.b

    def align(self):
        self.cost(len(self.A), len(self.B))
        self.RecurBackTrack(len(self.A), len(self.B))
        self.a = self.num_to_sequence(self.a)
        self.b = self.num_to_sequence(self.b)


### END OF CLASS BLOCK
################################################################################

print("Question 1: What is the optimal (here maximal) cost of an alignment of AATAAT\
 and AAGG using the above substitution matrix and gap cost -5?")
sequence1 = "AATAAT"
sequence2 = "AAGG"

myAlignment = Alignment([sequence1, sequence2])
myAlignment.align()
print(myAlignment.score)

print("Question 2: What is the optimal (here maximal) cost of an alignment of \
seq1.fasta and  seq2.fasta using the same substitution matrix and gap cost?")
sequence1 = "GGCCTAAAGGCGCCGGTCTTTCGTACCCCAAAATCTCGGCATTTTAAGATAAGTGAGTGTTGCGTTACACTAGCGATCTACCGCGTCTTATACTTAAGCGTATGCCCAGATCTGACTAATCGTGCCCCCGGATTAGACGGGCTTGATGGGAAAGAACAGCTCGTCTGTTTACGTATAAACAGAATCGCCTGGGTTCGCGGCCTAAAGGCGCCGGTCTTTCGTACCCCAAAATCTCGGCATTTTAAGATAAGTGAGTGTTGCGTTACACTAGCGATCTACCGCGTCTTATACTTAAGCGTATGCCCAGATCTGACTAATCGTGCCCCCGGATTAGACGGGCTTGATGGGAAAGAACAGCTCGTCTGTTTACGTATAAACAGAATCGCCTGGGTTCGC"
sequence2 = "GGGCTAAAGGTTAGGGTCTTTCACACTAAAGAGTGGTGCGTATCGTGGCTAATGTACCGCTTCTGGTATCGTGGCTTACGGCCAGACCTACAAGTACTAGACCTGAGAACTAATCTTGTCGAGCCTTCCATTGAGGGTAATGGGAGAGAACATCGAGTCAGAAGTTATTCTTGTTTACGTAGAATCGCCTGGGTCCGC"

myAlignment = Alignment([sequence1, sequence2])
myAlignment.align()
print(myAlignment.score)

print("Question 3: How does an optimal alignment look like for the above two \
pairs of sequences using the given substitution matrix and gap cost -5? ")
print(myAlignment.a)
print(myAlignment.b)


# Score matrix print
# for row in range(len(T)):
#     print(T[row])
