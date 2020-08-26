# PySeq: Python implementation of Needleman-Wunsch (N-W) algorithm and Hirschberg’s algorithm

Proteins consist of amino acid chains, and aligning the amino acid sequences of related proteins may show which regions have functional or structural importance. Global sequence alignment attempts to find the optimal alignment of two sequences of characters across their entire spans. In this project, we implement two dynamic programming algorithms for global sequence alignment: the Needleman-Wunsch algorithm and Hirschberg’s algorithm in Python. This is a course project that is intended to examine the theoretical runtime experimentally, through a sophisticated sequential implementation of the two algorithms in Python. 

**Development Team:**
* Mohammed Al Farhan (mohammed.farhan@kaust.edu.sa)
* Ahmad Shono (ahmad.shono@kaust.edu.sa)
* Liam Mencel (liam.mencel@kaust.edu.sa)
* Nabeelah Ali (nabeelah.ali@kaust.edu.sa)

## Algorithm Details
### Needleman-Wunsch (N-W) Algorithm

The N-W algorithm is a dynamic programming algorithm that builds up the best alignment using optimal alignments of smaller subsequences. This is achieved by filling all cells of a (n + 1, m + 1) matrix (where n and m are the lengths of the two sequences to be compared) according to the N-W recurrence relation and the chosen manipulation scores. 

**The algorithm can be implemented in three steps:**
*  Initialising the score matrix SM: The chosen substitution scores are stored in a score matrix SM, where each cell SM(xi, yi) reflects the score of the substitution of one character with another character from the specified character domain. The size of SM depends on the size of the character domain. For example, the size of SM is 44 for DNA with character domain C, T, A, G and is 2020 for proteins with character domain {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V}. For the purposes of this project, and for the simplicity, the BLOSUM62 (Appendix A) score matrix for proteins is used
*  Filling the traceback matrix TM and obtain T(n, m) as the similarity score: The traceback matrix TM is filled using the following recurrence relation: `TM(i,j) = max[TM(i-1, j) + gapPenalty, TM(i, j-1) + gapPenalty, TM(i − 1, j − 1) + SM(xj, yj )]`. Note that for the first row (i, j = 0) and the first column (i = 0, j), the formula will return a multiple of the gap penalty
*  Deducing the optimal alignment from the traceback matrix TM: Start from TM(n, m) and follow below condition rules until reaching TM(1, 1)
```
if TM(i, j) = TM(i − 1, j − 1) + SM(xi, yj )TM(i − 1, j − 1) then
  move diagonally (match/mismatch) and repeat on TM(i − 1, j − 1)
else if TM(i, j) = TM(i − 1, j) + gapPenalty then
  move left (gap in second sequence) and repeat on TM(i − 1, j)
else if TM(i, j) = TM(i, j − 1) + gapPenalty then
  move up (gap in first sequence) and repeat on TM(i, j − 1)
end if
```
### Hirschberg’s Algorithm

Hirschberg’s Algorithm can be described as a "divide and conquer" version of the Needleman-Wunsch algorithm. The key advantage of it is that it uses space complexity which is only linear in the lengths of the strings. In this algorithm, we will have a forwards subprogram and backwards subprogram, described below. In this method we will initialise the score matrix SM as with the N-W algorithm, which will be used to evaluate similarity between characters in the same way. Let score(x, y) denote the similarity score between two sequences x and y. Also, let pref ix(x, i) denote the subsequence of x which is comprised of the first i characters. Similarly, suff ix(x, i) denotes the subsequence comprised of the final i characters of x. Our algorithm can be implemented as a combination of three main routines:

**Forward(x, y)**

This routine takes as input two sequences x and y, and outputs an array of length m (size of y), which holds the m different scores between x and a prefix of y.
The routine starts by initialising an empty matrix T of size (n + 1)(m + 1) like in the N-W algorithm. Then, it sets T(0, j) = j*gapPenalty, for each 0 <= j <= m, where gapPenalty is the score applied for the insertion of a character. Now execute the following loop:
```
for i from 1 to n do
    T(i, 0) = T(i − 1, 0) + gapPenalty
    for j from 1 to m do
      T(i, j) = min[TM(i-1, j) + gapPenalty, TM(i, j-1) + gapPenalty, TM(i − 1, j − 1) + SM(xj, yj )]
      Delete T(i − 1, j − 1) from memory
    end for
  Delete T(i − 1, m) from memory
end for
```
Note that this loop caculates the matrix values row by row, similarly to the N-W algorithm. However, each time we store a new value we delete the corresponding value diagonally up-left from it, which is no longer needed. After the loop finishes we end up with the final row. Finally, it outputs row T(n, j), for 0 <= j <= m.

**Backward(x, y)**

This routine is essentially identical to the forwards routine, but will use suffixes of the strings instead of prefixes. It will again output an array of length m. The routine follows the following steps:

First, initialise an empty Matrix S of size (n + 1)(m + 1). This is almost the same as our previous matrix T but this time S(i, j) is the similarity scores between the suffix of x of size I and the suffix of y of size j.
Then, set S(0, j) = j*gapPenalty for each 0 <= j <= m, where gapPenalty is the scored applied for insertion of a character.
Now, execute the following loop:
```
for i from 1 to n do
    S(i, 0) = S(i − 1, 0) + gapPenalty
    for j from 1 to m do
      S(i, j) = min[S(i-1, j) + gapPenalty, S(i, j-1) + gapPenalty, S(i − 1, j − 1) + SM(n-i-1, n-j-1)]
      Delete S(i − 1, j − 1) from memory
    end for
  Delete S(i − 1, j − 1) from memory
end for
```
Note that this loop calculates the matrix values row by row, similarly to the N-W algorithm. However, each time we store a new value we delete the corresponding value diagonally up-left from it, which is no longer needed. After the loop finishes we end up with the final row.
Finally, Output row {S(n, j), for 0 <= j <= m}.

**Hirschberg(x, y)**

This recursive routine will call Forwards() and Backwards() as subroutines, and will accept as input the two sequences we wish to compute the score for. It will otuput the optimal alignment of x and y. 

If n <= 1 or m <= 1, we trivially output the alignment of x and y using the standard N-W algorithm. Else, partition x into two halves; *xlef* t is the first half of x, and *xright* is the final half. This is where the "Divide and Conquer" method plays its role. 
Next, we will call the Forwards and Backwards routines once each. Let F = Forwards(xlef t, y) and B = Backwards(right, y). So we have two arrays each of size m. 
Then, compute cut = value of j which minimises: F(j) + B(m − j); 0 <= j <= m. This can be computed in linear time; reading off the values from F and B for each value of j. Note that for each j, the quantity F(j) + B(m − j) represents the score of matching xlef t with pref ix(y, j), added to the score of matching xright with suff ix(y, m−j), i.e. the two partitioned sequences are matched in a way that forces the left half of x to be matched to the left partition of y, whilst the right half of x matches to the right partition of y. Note that the final optimal alignment will definitely divide this way for some value of j, hence, by finding the value of j which minimises the score we have in fact determined that the optimal alignment splits at cut. 
Now, we have a splitting point, we can recursively call the Hirschberg algorithm on the smaller subsequences of x and y, and output the corresponding alignments:

* Hirschberg(*xleft*, prefix(y, cut))
* Hirschberg(*xright*, suffix(y, m − cut))

By combining the alignment of the prefixes of x and y, with the alignment of the suffixes of x and y, we get an alignment between x and y. Finally, return concatenated alignment of x and y.

## Implementation Details

The scripts for this project were written in Python 2.7. Two functions were created, “nw” and “hirschberg”. Both take as input two string sequences, and output the optimal alignment and corresponding similarity score.

"nw" creates an (n + 1)(m + 1) matrix, implemented using the array data structure. This has the advantage of constant time storage and look-up of elements in the matrix. The general process of this function adheres to the recursive relation described in the N-W algorithm section.

In the "hirschberg" function, at each recursive step we call the "forwards" and "backwards" subroutines to determine the optimal splitting point in linear time. All data is stored in arrays, however it can be removed from memory once we find the splitting point. Once we determine the splitting point we divide the strings into two partitions each, and recall the "hirschberg" function on the pair of prefixes and pair of suffixes. When the recursion reaches the base cases, strings of length 0 or 1, we simply import the "nw" function described above and use it to solve the trivial cases.
