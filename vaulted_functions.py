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

    return a, b
