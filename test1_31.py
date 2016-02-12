import random

# import our algorithms for testing
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

def alignIndex(alignment, i):
    # returns index of i+1th non-gap character in string
    pointer = 0
    charsPassed = 0
    while charsPassed < i+1:
        if not (alignment[pointer] == '-'):
            charsPassed += 1
        pointer += 1
    return pointer-1

def mutate(string, steps):
    # Will mutate a string using a fixed number of edits
    parent = parentAlign = child = childAlign = string
    for r in range(steps):
        m = random.randint(1, 3)
        if not child:
            m = 2
        if (m==1):
            #Delete
            i = random.randint(0, len(child)-1)
            j = alignIndex(childAlign, i)
            child = child[:i] + child[i+1:]
            childAlign = childAlign[:j] + '-' + childAlign[j+1:]
        elif (m==2):
            #Insert
            i = random.randint(0, len(child))
            char = random.choice(alphabet)
            if(i == len(child)):
                child = child + char
                childAlign = childAlign + char
                parentAlign = parentAlign + '-'
            else:
                j = alignIndex(childAlign, i)
                child = child[:i] + char + child[i:]
                childAlign = childAlign[:j] + char + childAlign[j:]
                parentAlign = parentAlign[:j] + '-' + parentAlign[j:]
        elif (m==3):
            #Substitute
            i = random.randint(0, len(child)-1)
            char = random.choice(alphabet)
            j = alignIndex(childAlign, i)
            child = child[:i] + char + child[i+1:]
            childAlign = childAlign[:j] + char + childAlign[j+1:]
    # Output mutated string, and score applied
    score = 0
    for i in range(len(parentAlign)):
        if(parentAlign[i] == '-' and childAlign[i] == '-'):
            continue
        elif(parentAlign[i] == '-' or childAlign[i] == '-'):
            score += gapPenalty
        else:
            score += simMatrix[alphEnum[parentAlign[i]]][alphEnum[childAlign[i]]]
    return (child, score)  
        
# Change this seed to control deterministic generation of random sequences
random.seed()

# Number of trials
trials = 1000

print "Testing..."

for t in range(trials):
    fail = False
    n = random.randint(5, 100)
    k = random.randint(1, n)
    parent = ""
    for i in range(n):
        parent = parent + random.choice(alphabet)
    child, work = mutate(parent, k)
    nwout = nw(parent, child, simMatrix, gapPenalty, alphEnum)
    hbout = hirschberg(parent, child, simMatrix, gapPenalty, alphEnum)
    
    #Verification tests
    compress = nwout[0].replace('-', '')
    if compress != parent:
        fail = True
    compress = nwout[1].replace('-', '')
    if compress != child:
        fail = True
    compress = hbout[0].replace('-', '')
    if compress != parent:
        fail = True
    compress = hbout[1].replace('-', '')
    if compress != child:
        fail = True
        
    score = 0
    for i in range(len(nwout[0])):
        if nwout[0][i] == '-' and nwout[1][i] == '-':
            fail = True
            break
        elif nwout[0][i] == '-':
            score += gapPenalty
        elif nwout[1][i] == '-':
            score += gapPenalty
        else:
            score += simMatrix[alphEnum[nwout[0][i]]][alphEnum[nwout[1][i]]]
    if score != nwout[2]:
        fail = True
    score = 0
    for i in range(len(hbout[0])):
        if hbout[0][i] == '-' and hbout[1][i] == '-':
            fail = True
            break
        elif hbout[0][i] == '-':
            score += gapPenalty
        elif hbout[1][i] == '-':
            score += gapPenalty
        else:
            score += simMatrix[alphEnum[hbout[0][i]]][alphEnum[hbout[1][i]]]
    if score != hbout[2]:
        fail = True

    if fail:
        g = open("errorlist.txt", 'w')
        g.write("FAIL - Verification\n")
        g.write(parent + '\n')
        g.write(child + "\n\n")
        g.close()
        print "Verification error - see text file"
        break

    #Bounds test
    if nwout[2] != hbout[2] or nwout[2] < work or hbout < work:
        g = open("errorlist.txt", 'w')
        g.write("FAIL - Bounds\n")
        g.write(parent + '\n')
        g.write(child + "\n\n")
        g.close()
        print "Bounds error - see text file"
        break
else:
    g = open("errorlist.txt", 'w')
    g.write("No errors detected")
    g.close()
    print "Done, no errors detected"
