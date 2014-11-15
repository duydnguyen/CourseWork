#!/user/bin/env python
# Author: Duy Nguyen
# Solution to question 3,  homework 4
# Viterbi Algorithm
# You should be able to run the program as of the following commands


import sys
import io
import getopt


def read_transitions(filename):
    "Read transition probabilities between states from transitions.txt"
    ### NEED TO FIX the case when non-consecutive states
    probs = []
    states = []
    i = 0
    ii = 0
    with open(filename, 'r') as file:
        for line in file: 
            #import pdb; pdb.set_trace()
            vals = line.strip().split()
            states.insert(i, int(vals[0]) )
            states.insert(i+1, int(vals[1]) )
            probs.insert(ii, float(vals[2]))
            i += 2
            ii += 1
    numStates = len(set(states))        
    ## Create  the transitions probability matrix filled with -1
    trans = [[-1 for x in range(numStates)] for y in range(numStates)]
    # fill the matrix trans
    ii = 0
    for i in range(len(states)):
        if (i % 2 == 0):
            row = states[i]
            col = states[i+1]
            trans[row][col] = probs[ii]
            ii += 1
    return trans

def read_emissions(filename, numStates):
    "Read emission probabilities from emissions.txt with order A, C, G, Tx"
    ### NEED to fix the case when states are not in consecutive orders; capital/noncapital A, C, G, T
    ### map state 0 -> begin, state numState -> end
    emissions = {}
    index = {'A':1, 'C':2, 'G':3, 'T':4}
    # Initilize emission dictionary
    for i in range(numStates-2):
        emissions[i+1] = [-1,-1,-1,-1]
    with open(filename, 'r') as file:
        for line in file:
            vals = line.strip().split()
            state = int(vals[0])
            letter = vals[1]
            pos = index[letter]
            prob = float(vals[2])
            # Add to dictionary emissions
            emissions[state][pos-1] = prob 
    return emissions

def multList(L1, L2):
    return [x*y for x, y in zip(L1, L2)]

def Initilize_V(numStates, seqLen):
    "Initilize matrix v for Viterbi algorithm: row = states, column = t (seq. position) "
    v = []
    v = [[-1 for x in range(seqLen+1)] for y in range(numStates)]
    # let row 0 (begin state) = 0's
    v[0] = [0 for x in range(seqLen+1)]
    v[0][0] = 1
    # let row = numStates (end state) = 0's
    for i in range(0, seqLen+1):
        v[numStates-1][i] = 0    
    # let col 0 (t = 0) = 0's
    for i in range(1, numStates):
        v[i][0] = 0
    return v

def transpose_List(l):
    return  [list(i) for i in zip(*l)]

def evalViterbi(transitions, emission, begin_state, end_state, sequence):
    "Viterbi algorithm"
    ### NEED to handle case sequence not capitalized
    index = {'A':1, 'C':2, 'G':3, 'T':4}
    seqLen = len(sequence)
    numStates = len(transitions)
    # Transpose of transitions matrix
    trans_t = transpose_List(transitions)
    # Initialize matrix v
    v = Initilize_V(numStates, seqLen)
    ## Main 
#    import pdb; pdb.set_trace()
    for t in range(1, seqLen+1):
        x_t = index[sequence[t-1]]
        # only handle case when states are not begin state; does not compute end state 
        for k in range(1, numStates-1):
            trans_v = transpose_List(v)
            evalMax = max(multList(trans_t[k], trans_v[t-1]))
            emiss = emissions[k][x_t-1]
            v[k][t] = emiss * evalMax
    print(v)
    return 0 
    
if __name__ == '__main__':
    transitions = read_transitions('Tests/transition2.txt')
    numStates = len(transitions)
    emissions = read_emissions('Tests/emission2.txt', numStates)
    begin_state = 0
    end_state = 5
    sequence = 'TAG'
    evalViterbi(transitions, emissions, begin_state, end_state, sequence)
    

    
