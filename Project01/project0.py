#Implementing alignment with dynamic programming

def sequence_to_num(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', '0')
    sequence = sequence.replace('C', '1')
    sequence = sequence.replace('G', '2')
    sequence = sequence.replace('T', '3')

    return sequence

def num_to_sequence(number):
    number = number.upper()
    number = number.replace('0', 'A')
    number = number.replace('1', 'C')
    number = number.replace('2', 'G')
    number = number.replace('3', 'T')

    return number

def cost(i,j,k=-5):

    if T[i][j] == float("-inf"):
        v1 = v2 = v3 = v4 = float("-inf")

        if (i > 0) and (j > 0):
            v1 = cost(i-1, j-1, k) + N[int(A[i-1])][int(B[j-1])]
        if (i > 0) and (j >= 0):
            v2 = cost(i-1, j, k) + k
        if (i >= 0) and (j > 0):
            v3 = cost(i, j-1, k) + k
        if (i == 0) and (j == 0):
            v4 = 0

        T[i][j] = max(v1,v2,v3,v4)

    return T[i][j]

def RecurBackTrack(i, j, a, b, k=-5):

    if(i > 0) and (j > 0) and T[i][j] == (T[i-1][j-1] + N[int(A[i-1])][int(B[j-1])]):
        a = A[i-1] + a
        b = B[j-1] + b
        return RecurBackTrack(i-1, j-1, a, b)
    elif(i > 0) and (j >= 0) and T[i][j] == T[i-1][j] + k:
        a = A[i-1] + a
        b = "-" + b
        return RecurBackTrack(i-1, j, a, b)
    elif(i >= 0) and (j > 0) and T[i][j] == T[i][j-1] + k:
        a = "-" + a
        b = B[j-1] + b
        return RecurBackTrack(i, j-1, a, b)

    return a, b

### END OF FUNCTIONS BLOCK
################################################################################

N = [[10,2,5,2],
    [2,10,2,5],
    [5,2,10,2],
    [2,5,2,10]] # substitution maxtrix (given)
k = -5      # gap cost (given)
# sequence1 = "GGCCTAAAGGCGCCGGTCTTTCGTACCCCAAAATCTCGGCATTTTAAGATAAGTGAGTGTTGCGTTACACTAGCGATCTACCGCGTCTTATACTTAAGCGTATGCCCAGATCTGACTAATCGTGCCCCCGGATTAGACGGGCTTGATGGGAAAGAACAGCTCGTCTGTTTACGTATAAACAGAATCGCCTGGGTTCGCGGCCTAAAGGCGCCGGTCTTTCGTACCCCAAAATCTCGGCATTTTAAGATAAGTGAGTGTTGCGTTACACTAGCGATCTACCGCGTCTTATACTTAAGCGTATGCCCAGATCTGACTAATCGTGCCCCCGGATTAGACGGGCTTGATGGGAAAGAACAGCTCGTCTGTTTACGTATAAACAGAATCGCCTGGGTTCGC"
# sequence2 = "GGGCTAAAGGTTAGGGTCTTTCACACTAAAGAGTGGTGCGTATCGTGGCTAATGTACCGCTTCTGGTATCGTGGCTTACGGCCAGACCTACAAGTACTAGACCTGAGAACTAATCTTGTCGAGCCTTCCATTGAGGGTAATGGGAGAGAACATCGAGTCAGAAGTTATTCTTGTTTACGTAGAATCGCCTGGGTCCGC"
sequence1 = "CGTGTCAAGTCT"
sequence2 = "ACGTCGTAGCTAGG"
T = [ [ float("-inf") for i in range(len(B)+1) ] for j in range(len(A)+1) ] #inizilize score matrix

# A and B are sequences translated to numbers (0-3)
A,B = sequence_to_num(sequence1), sequence_to_num(sequence2)

res = cost(len(A), len(B), k)
res

# Score matrix print
for row in range(len(T)):
    print(T[row])

a = ""
b = ""
alignment_numbers = RecurBackTrack(len(A), len(B), a, b, k)
alignment = [num_to_sequence(alignment_numbers[0]), num_to_sequence(alignment_numbers[1])]
print(alignment[0])
print(alignment[1])

### TRY HERE
################################################################################
