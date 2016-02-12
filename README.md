# Sequence Alignment: A Python Implementation of Needleman-Wunsch (N-W) Algorithm as well as Hirschberg’s Algorithm

Proteins consist of amino acid chains, and aligning the amino acid sequences of related proteins may show which regions have functional or structural importance. Global sequence alignment attempts to find the optimal alignment of two sequences of characters across their entire spans. In this project, we implement two dynamic programming algorithms for global sequence alignment: the Needleman-Wunsch algorithm and Hirschberg’s algorithm in Python. This is a course project that is intended to examine the theoretical runtime experimentally, through a sophisticated sequential implementation of the two algorithms in Python.

**Development Team:**
* Liam Mencel (liam.mencel@kaust.edu.sa)
* Mohammed Al Farhan (mohammed.farhan@kaust.edu.sa)
* Ahmad Shono (ahmad.shono@kaust.edu.sa)
* Nabeelah Ali (nabeelah.ali@kaust.edu.sa)

## Algorithm Details
### Needleman-Wunsch (N-W) Algorithm

The N-W algorithm is a dynamic programming algorithm that builds up the best alignment using optimal alignments of smaller subsequences. This is achieved by filling all cells of a (n + 1, m + 1) matrix (where n and m are the lengths of the two sequences to be compared) according to the N-W recurrence relation and the chosen manipulation scores. 

**The algorithm can be implemented in three steps:**
*  Initialising the score matrix SM: The chosen substitution scores are stored in a score matrix SM, where each cell SM(xi, yi) reflects the score of the substitution of one character with another character from the specified character domain. The size of SM depends on the size of the character domain. For example, the size of SM is 44 for DNA with character domain C, T, A, G and is 2020 for proteins with character domain {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V}. For the purposes of this project, and for the simplicity, the BLOSUM62 (Appendix A) score matrix for proteins is used
*  Filling the traceback matrix TM and obtain T(n, m) as the similarity score: The traceback matrix TM is filled using the following recurrence relation: `TM(i,j) = max[TM(i-1, j) + gapPenalty, TM(i, j-1) + gapPenalty, TM(i − 1, j − 1) + SM(xj, yj )]` Note that for the first row (i, j = 0) and the first column (i = 0, j), the formula will return a multiple of the gap penalty
*  Deducing the optimal alignment from the traceback matrix TM: Start from TM(n, m) and follow below condition rules until reaching TM(1, 1)

`if TM(i, j) = TM(i − 1, j − 1) + SM(xi, yj )TM(i − 1, j − 1) then`
` move diagonally (match/mismatch) and repeat on TM(i − 1, j − 1)`
`else if TM(i, j) = TM(i − 1, j) + gapP enalty then`
` move left (gap in second sequence) and repeat on TM(i − 1, j)`
`else if TM(i, j) = TM(i, j − 1) + gapP enalty then`
` move up (gap in first sequence) and repeat on TM(i, j − 1)`
`end if`
