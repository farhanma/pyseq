from functions1_3 import nw
from functions1_3 import forwards
from functions1_3 import backwards
from functions1_3 import hirschberg

# Read alphabet and scores from text file
f = open("scores.txt", 'r')
alphabet = f.readline()
while(alphabet[-1] in ['\n', '\r']):
    alphabet = alphabet[:-1]
f.readline()
gapPenalty = int(f.readline())
f.readline()
simMatrix = []
line = f.readline()
while(line):
    row = list(int(x) for x in line.split())
    simMatrix.append(row)
    line = f.readline()
f.close()

# Create a 1-1 mapping from characters to integers, for simplicity in algorithm
alphEnum = dict([(alphabet[i], i) for i in range(len(alphabet))])

# Load input sequences
f = open("sequences.txt", 'r')
line = f.readline()

# Open output file, in preparation for storing output alignments
g = open("alignments.txt", 'w')

while(line):
    # This loop repeats until no more input sequences are found
    # At each iteration we read the next two sequences and run the algorithm on them
    A = line
    while(A[-1] in ['\n', '\r']):
        A = A[:-1]
    B = f.readline()
    while(B[-1] in ['\n', '\r']):
        B = B[:-1]
    f.readline()
    line = f.readline()

    # Run the Hirschberg algorithm

    print "First sequence:", A
    print "Second sequence:", B
    print "Calculating alignment distance by Hirschberg method..."
    
    z = hirschberg(A, B, simMatrix, gapPenalty, alphEnum)

    print "Alignment of A: ", z[0]
    print "Alignment of B: ", z[1]
    print "Similarity score: ", z[2], '\n'

    # Write outputs to text file
    g.write(str(z[2]) + "\n")
    g.write(z[0] + "\n")
    g.write(z[1] + "\n")
    g.write("\n")

# Close the files and finish
f.close()
g.close()
