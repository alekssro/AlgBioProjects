def cost(self, i, j):
    """Recursive function to perform a constant gap cost alignment"""

    k = self.k
    if self.T[i][j] == float("inf"):
        v1 = v2 = v3 = v4 = float("inf")

        if (i > 0) and (j > 0):
            v1 = self.cost(i-1, j-1) - self.N[int(self.A[i-1])][int(self.B[j-1])]
        if (i > 0) and (j >= 0):
            v2 = self.cost(i - 1, j) + self.gapcost(k)
        if (i >= 0) and (j > 0):
            v3 = self.cost(i, j - 1) + self.gapcost(k)
        if (i == 0) and (j == 0):
            v4 = 0

        self.T[i][j] = min(v1, v2, v3, v4)

    self.score = self.T[i][j]
    return self.T[i][j]

def cost_iter(self):
    """Iterative function to perform a constant gap cost alignment"""

    # First column (0)
    self.T[0][0] = 0
    for i in range(1, len(self.T)):
        self.T[i][0] = self.T[i-1][0] + self.gapcost(i)
        # print(self.T[i][0])

    # Rest of the matrix
    for j in range(1, len(self.T[0])):
        v1 = v2 = v3 = float("inf")
        for i in range(0, len(self.T)):

            if (i > 0) and (j > 0):     # Diagonal arrow
                v1 = self.T[i-1][j-1] + self.N[int(self.A[i-1])][int(self.B[j-1])]
            if (i > 0) and (j >= 0):    # Vertical arrow
                v2 = self.T[i - 1][j] + self.gapcost(j)
            if (i >= 0) and (j > 0):    # Horizontal arrow
                v3 = self.T[i][j - 1] + self.gapcost(i)

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

def affine_align(self):
    """Global alignment with affine or linear penalties. We assume we are minimizing."""
    self.D = self.initMatrix(len(self.seq1) + 1, len(self.seq2) + 1, value=0)
    self.Ins = self.initMatrix(len(self.seq1) + 1, len(self.seq2) + 1, value=0)

    self.T[0][0] = 0
    # Rest of the matrix
    for j in range(0, len(self.T[0])):
        for i in range(1, len(self.T)):
            v1 = v2 = v3 = float("inf")

            if (i > 0) and (j > 0):     # Diagonal arrow
                v1 = self.T[i-1][j-1] + self.N[int(self.A[i-1])][int(self.B[j-1])]
            if (i > 0) and (j >= 0):    # Vertical arrow
                d1 = d2 = float("inf")
                d1 = self.T[i-1][j] + self.gapcost(1)
                if (i > 1) and (j >= 0):
                    d2 = self.D[i-1][j] + self.gapcost(0)
                self.D[i][j] = min(d1, d2)
                v2 = self.D[i][j]
            if (i >= 0) and (j > 0):    # Horizontal arrow
                i1 = i2 = float("inf")
                i1 = self.T[i][j-1] + self.gapcost(1)
                if (i >= 0) and (j > 1):
                    i2 = self.Ins[i][j-1] + self.gapcost(0)
                self.Ins[i][j] = min(i1, i2)
                v3 = self.Ins[i][j]

            self.T[i][j] = min(v1, v2, v3)

    self.score = self.T[i][j]

    return self.T[i][j]
