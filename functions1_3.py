def nw(A, B, simMatrix, gapPenalty, alphEnum):
    # The Needleman-Wunsch algorithm
    
    # Stage 1: Create a zero matrix and fills it via algorithm
    n, m = len(A), len(B)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(m+1):
        mat[0][j] = gapPenalty*j
    for i in range(n+1):
        mat[i][0] = gapPenalty*i
    for i in range(1, n+1):
        for j in range(1, m+1):
            mat[i][j] = max(mat[i-1][j-1] + simMatrix[alphEnum[A[i-1]]][alphEnum[B[j-1]]], mat[i][j-1] + gapPenalty, mat[i-1][j] + gapPenalty)

    # Stage 2: Computes the final alignment, by backtracking through matrix
    alignmentA = ""
    alignmentB = ""
    i, j = n, m
    while i and j:
        score, scoreDiag, scoreUp, scoreLeft = mat[i][j], mat[i-1][j-1], mat[i-1][j], mat[i][j-1]
        if score == scoreDiag + simMatrix[alphEnum[A[i-1]]][alphEnum[B[j-1]]]:
            alignmentA = A[i-1] + alignmentA
            alignmentB = B[j-1] + alignmentB
            i -= 1
            j -= 1
        elif score == scoreUp + gapPenalty:
            alignmentA = A[i-1] + alignmentA
            alignmentB = '-' + alignmentB
            i -= 1
        elif score == scoreLeft + gapPenalty:
            alignmentA = '-' + alignmentA
            alignmentB = B[j-1] + alignmentB
            j -= 1
    while i:
        alignmentA = A[i-1] + alignmentA
        alignmentB = '-' + alignmentB
        i -= 1
    while j:
        alignmentA = '-' + alignmentA
        alignmentB = B[j-1] + alignmentB
        j -= 1
    # Now return result in format: [1st alignment, 2nd alignment, similarity]
    return [alignmentA, alignmentB, mat[n][m]]

def forwards(x, y, simMatrix, gapPenalty, alphEnum):
    # This is the forwards subroutine.
    n, m = len(x), len(y)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(m+1):
        mat[0][j] = gapPenalty*j
    for i in range(1, n+1):
        mat[i][0] = mat[i-1][0] + gapPenalty
        for j in range(1, m+1):
            mat[i][j] = max(mat[i-1][j-1] + simMatrix[alphEnum[x[i-1]]][alphEnum[y[j-1]]],
                            mat[i-1][j] + gapPenalty,
                            mat[i][j-1] + gapPenalty)
        # Now clear row from memory.
        mat[i-1] = []
    return mat[n]    

def backwards(x, y, simMatrix, gapPenalty, alphEnum):
    # This is the backwards subroutine.
    n, m = len(x), len(y)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(m+1):
        mat[0][j] = gapPenalty*j
    for i in range(1, n+1):
        mat[i][0] = mat[i-1][0] + gapPenalty
        for j in range(1, m+1):
            mat[i][j] = max(mat[i-1][j-1] + simMatrix[alphEnum[x[n-i]]][alphEnum[y[m-j]]],
                            mat[i-1][j] + gapPenalty,
                            mat[i][j-1] + gapPenalty)
        # Now clear row from memory.
        mat[i-1] = []
    return mat[n]

def hirschberg(x, y, simMatrix, gapPenalty, alphEnum):
    # This is the main Hirschberg routine.
    n, m = len(x), len(y)
    if n<2 or m<2:
        # In this case we just use the N-W algorithm.
        return nw(x, y, simMatrix, gapPenalty, alphEnum)
    else:
        # Make partitions, call subroutines.
        F, B = forwards(x[:n/2], y, simMatrix, gapPenalty, alphEnum), backwards(x[n/2:], y, simMatrix, gapPenalty, alphEnum)
        partition = [F[j] + B[m-j] for j in range(m+1)]
        cut = partition.index(max(partition))
        # Clear all memory now, so that we don't store data during recursive calls.
        F, B, partition = [], [], []
        # Now make recursive calls.
        callLeft = hirschberg(x[:n/2], y[:cut], simMatrix, gapPenalty, alphEnum)
        callRight = hirschberg(x[n/2:], y[cut:], simMatrix, gapPenalty, alphEnum)
        # Now return result in format: [1st alignment, 2nd alignment, similarity]
        return [callLeft[r] + callRight[r] for r in range(3)]
