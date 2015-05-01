#!/user/bin/env python
# Author: Duy Nguyen
# BMI/CS_776: Nussinov Algorithm
# You should be able to run the program as of the following command
# python Nussinov.py RNASEQUENCE

import sys
import io
import getopt
import numpy as np
# import doctest
from math import log

def read_sequences(filename):
    " Read sequences from file "
    sequences = []
    with open(filename, 'r') as file:
        for line in file:
            seq = line.strip().split()
            sequences.append(seq)
    return sequences

def evaldelta(seq, i, j):
    """Evaluate delta function: 1(x_i and x_j are complementary); i ,j in 1:L
    
    >>> evaldelta('UCCAGG', 1, 4)
    1
    >>> evaldelta('UCCAGG', 1, 2)
    0
    >>> evaldelta('UCCAGG', 3, 6)
    1
    """
    delta = 0
    dict = {'A':'U', 'C':'G', 'G':'C', 'U':'A'}
    x_i = seq[i-1]
    x_j = seq[j-1]
    if (dict[x_i] == x_j):
        delta = 1
    return delta

def display(twolist):
    "Display a 2-dim list as an 2-dim array"
    L = len(twolist)
    for i in range(L):
        print(twolist[i])
    return 0

def createPath(L):
    """Create the path which starts at the lower right conner of  Gamma matrix
    
    >>> createPath(6)
    [(5, 6), (4, 5), (4, 6), (3, 4), (3, 5), (3, 6), (2, 3), (2, 4), (2, 5), (2, 6), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6)]
    """
    path = []
    for i in range(1,L)[::-1]:
        for j in range(i+1, L+1):
            path.append((i,j))
    return path

def evalGamma(seq, L):
    """Eval matrix Gamma for Nussinov Algorithm
    
    >>> evalGamma('UCCAGG', 6)
    [['*', '*', '*', '*', '*', '*', '*'], ['*', 0, 0, 0, 1, 1, 2], ['*', 0, 0, 0, 0, 1, 2], ['*', '*', 0, 0, 0, 1, 1], ['*', '*', '*', 0, 0, 0, 0], ['*', '*', '*', '*', 0, 0, 0], ['*', '*', '*', '*', '*', 0, 0]]
    >>> evalGamma('AAAUCCCAGGA', 11)
    [['*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*'], ['*', 0, 0, 0, 1, 1, 1, 1, 1, 2, 3, 3], ['*', 0, 0, 0, 1, 1, 1, 1, 1, 2, 3, 3], ['*', '*', 0, 0, 1, 1, 1, 1, 1, 2, 3, 3], ['*', '*', '*', 0, 0, 0, 0, 0, 1, 1, 2, 3], ['*', '*', '*', '*', 0, 0, 0, 0, 0, 1, 2, 2], ['*', '*', '*', '*', '*', 0, 0, 0, 0, 1, 2, 2], ['*', '*', '*', '*', '*', '*', 0, 0, 0, 1, 1, 1], ['*', '*', '*', '*', '*', '*', '*', 0, 0, 0, 0, 0], ['*', '*', '*', '*', '*', '*', '*', '*', 0, 0, 0, 0], ['*', '*', '*', '*', '*', '*', '*', '*', '*', 0, 0, 0], ['*', '*', '*', '*', '*', '*', '*', '*', '*', '*', 0, 0]]
    """
    Gamma = [['*' for x in range(L+1)] for y in range(L+1)] 
    path = []
    delta = []
    path = createPath(L)
    # Initialization
    for i in range(2, L+1):
        Gamma[i][i-1] = 0
    for i in range(1, L+1):
        Gamma[i][i] = 0
    # Recursion
    # import pdb; pdb.set_trace()
    for entry in path:
        all_cases = []
        i = entry[0]
        j = entry[1]
        delta = evaldelta(seq, i, j)
        # print i, j, delta
        bifurcation = 0
        # Evaluate bifurcation case
        if j > i+1:
            # print i, j
            bif_list = []
            for k in range(i+1, j):
                bif_list.append(Gamma[i][k] + Gamma[k+1][j] )
                # print bif_list
            bifurcation = max(bif_list)
        all_cases = [Gamma[i+1][j], Gamma[i][j-1], Gamma[i+1][j-1] + delta, bifurcation ]
        # print all_cases
        Gamma[i][j] = max(all_cases)
    return Gamma


def traceback(seq ,Gamma):
    "Traceback solution"
    solution = []
    L = len(Gamma) - 1
    stack = [(1,L)]
    while stack:
        # import pdb; pdb.set_trace()
        delta  = 0
        #print stack
        entry = stack.pop()
        #print entry
        i = entry[0]
        j = entry[1]
        delta = evaldelta(seq, i, j)
        if j <= i:
            continue
        elif Gamma[i+1][j] == Gamma[i][j]:
            stack.append((i+1,j))
        elif Gamma[i][j-1] == Gamma[i][j]:
            stack.append((i,j-1))
        elif Gamma[i+1][j-1] + delta == Gamma[i][j]:
            solution.append((i,j))
            stack.append((i+1,j-1))
        else:
            for k in range(i+1, j):
                if Gamma[i][k] + Gamma[k+1][j] == Gamma[i][j]:
                    stack.append((k+1,j))
                    stack.append((i,k))
                    break
    return solution

def main(argv):
    """ Take arguments from the following command
    python Nussinov.py RNASEQUENCE
    """
    try:
        opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    sequence = args[0]
    return sequence

def output(solution):
    "Output traceback in lexicographical order"
    solution.sort()
    for entry in solution:
        print entry[0], entry[1]
    return 0
    

if __name__ == '__main__':
    # Input
    sequence = main(sys.argv[1:])
    # Nussinov Algorithm
    Gamma = []
    solution = []
    L = len(sequence)
    Gamma = evalGamma(sequence, L)
    solution = traceback(sequence, Gamma)
    output(solution)
    # doctest.testmod()

    

    
