#!/user/bin/env python
# Author: Duy Nguyen
# BMI/CS_776, Homework 3, Question 2
# PWM Models
# You should be able to run the program as of the following command
# python ???

import sys
import io
import getopt
import numpy as np
import doctest
from math import log

def read_sequences(filename):
    " Read sequences from file "
    sequences = []
    with open(filename, 'r') as file:
        for line in file:
            seq = line.strip().split()
            sequences.append(seq)
    return sequences

def normalize(list):
    " Normalize a 1-dim list"
    tot = sum(list)
    for i in range(len(list)):
        list[i] = float(list[i]) / tot
    return list

def get_column(mat, i):
    " get the ith column of a given matrix"
    return [row[i] for row in mat]

def count_Char(letter, sequence):
    """ Count the number of occurences of the given letter in sequence 

    >>> count_Char('A', 'GTCATTACTA')
    3
    >>> count_Char('C', 'GTCATTACTA')
    2
    >>> count_Char('G', 'GTCATTACTA')
    1
    >>> count_Char('T', 'GTCATTACTA')
    4
    """
    count = 0
    for i in range(len(sequence)):
         if (sequence[i] == letter):
             count += 1
    return count

def extract_col(sequences, col_index):
    """Given a set of sequences and column index, extract the corresponding column (in string)
    
    >>> extract_col([['GAGGTAAAC'],['TCCGTAAGT'],['CAGGTTGGA'],['ACAGTCAGT'],['TAGGTCATT'],['TAGGTACTG'],['ATGGTAACT'],['CAGGTATAC'],['TGTGTGAGT'],['AAGGTAAGT']], 0)
    'GTCATTACTA'
    """
    num_seqs = len(sequences)
    col_seq = ''
    for i in range(num_seqs):
        seq = sequences[i][0]
        #col_seq.append(seq[col_index])
        col_seq += seq[col_index]
    return col_seq
    
def learnPFM(sequences):
    """ Learn PFM (position frequency matrix) from sequence examples
    
    >>> learnPFM([['GAGGTAAAC'],['TCCGTAAGT'],['CAGGTTGGA'],['ACAGTCAGT'],['TAGGTCATT'],['TAGGTACTG'],['ATGGTAACT'],['CAGGTATAC'],['TGTGTGAGT'],['AAGGTAAGT']])
    [[3, 6, 1, 0, 0, 6, 7, 2, 1], [2, 2, 1, 0, 0, 2, 1, 1, 2], [1, 1, 7, 10, 0, 1, 1, 5, 1], [4, 1, 1, 0, 10, 1, 1, 2, 6]]
    """
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    ncol = len(sequences[0][0])
    mat =  [[0 for x in range(ncol)] for y in range(4)] 
    #import pdb; pdb.set_trace()
    for i in range(ncol):
        seq = extract_col(sequences, i)
        for j in range(4):
            mat[j][i] = count_Char(index[j], seq) 
    return mat

def learnPWM(PFM):
    """  Given a position frequency matrix, the function evaluates PWM (position weight matrix)
    """
    ncol = len(PFM[0])
    mat =  [[0 for x in range(ncol)] for y in range(4)]
    for col in range(ncol):
        seq_col = get_column(PFM, col)
        # plus 1 for pseudocount to avoid zeros
        for i in range(4):
            seq_col[i] += 1
        seq_col_normalized = normalize(seq_col)
        for i in range(4):
            mat[i][col] = seq_col_normalized[i]
    return mat


def evalScore(sequences, PWM_positive, PWM_negative):
    " Evaluate score(sequence) = log(P(seq|model_p)/P(seq|model_n))"
    scores = []
    numSeq = len(sequences)
    ncol = len(sequences[0][0])
    for i in range(numSeq):
        #import pdb; pdb.set_trace()
        seq = sequences[i][0]
        prob_positive = eval_prob_Zp(seq, 1, ncol, ncol, PWM_positive)
        prob_negative = eval_prob_Zp(seq, 1, ncol, ncol, PWM_negative)
        scores.append(log(prob_positive) - log(prob_negative))
    return scores

def find_Indices(sequence, char):
    """Given a sequence and a character in {A, C, G, T}, output the position of char in the seq
    >>> find_Indices('TTAAAATA', 'A')
    [2, 3, 4, 5, 7]

    >>> find_Indices('TTGTTGCC', 'T')
    [0, 1, 3, 4]
    """
    pos = []
    index = sequence.find(char)
    while (index != -1):
        pos.append(index)
        index = sequence.find(char, index + 1)
    return pos

def eval_Matches(col_i, col_j, char):
    """ Given a char in {A,C,G,T}, evaluate matches/non-matches between two sequence columns
    >>> eval_Matches('GTTTGCTA', 'TTAAAATA', 'A')
    [2, 3]

    >>> eval_Matches('GTTTGCTA', 'TTAAAATA', 'T')
    [2, 1]

    >>> eval_Matches('GTTTGCTA', 'TTGTTGCC', 'T')
    [4, 0]

    >>> eval_Matches('GTTTGCTA', 'TTGTTGCC', 'C')
    [0, 2]
    """
    pos = find_Indices(col_j, char)
    combined = []
    for i in range(len(pos)):
        s = col_i[pos[i]] + col_j[pos[i]]
        combined.append(s)
    # eval the number of matches
    matches = 0
    check = [True for i in range(len(combined))]
    for k in range(len(combined)):
        match_indices = [i for i, j in enumerate(combined) if j == combined[k] ]
        if (len(match_indices) > 1):
            for i in match_indices:
                if check[i]:
                    matches += 1
                    check[i] = False
    non_matches = len(pos) - matches
    return [matches, non_matches]

def eval_ChiSq(table):
    """ Given a 4-by-2 contigency table, evaluate the Chi Square Statistics
    >>> eval_ChiSq([[4, 2], [0, 2], [0, 1], [0, 1]])
    4.444444444444445
    """
    Chi = 0
    table_expected = [[0 for x in range(2)] for y in range(4)] 
    Row = []
    Col = []
    ## Compute the expected number of counts E_ij
    for i in range(4):
        Row.append(sum(table[i]))
    for j in range(2):
        Col.append(sum(get_column(table, j)))
    N = sum(Row)
    for i in range(4):
        for j in range(2):
            table_expected[i][j] = float(Row[i] * Col[j]) / N
            Chi += float(table[i][j] - table_expected[i][j] ) ** 2 / table_expected[i][j]
    return Chi
def eval_Table(sequences, index_i, index_j):
    """ Given set of sequences and indices i, j , eval a 4-by-2 contigency table
    >>> eval_Table( [['GAGGTAAAC'],['TCCGTAAGT'],['CAGGTTGGA'],['ACAGTCAGT'],['TAGGTCATT'],['TAGGTACTG'],['ATGGTAACT'],['CAGGTATAC'],['TGTGTGAGT'],['AAGGTAAGT']], 0, 1)
    [[4, 2], [0, 2], [0, 1], [0, 1]]
    """
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    table =  [[0 for x in range(2)] for y in range(4)]    
    col_i = extract_col(sequences, index_i )
    col_j = extract_col(sequences, index_j )

    for i in range(4):
        char = index[i]
        table[i] = eval_Matches(col_i, col_j, char)

    return table

def find_MDD_subtree(T, P):
    " Find the tree using the Maximal Dependence Decomposition (MDD) algorithm"
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    for i in P:
        PFM = learnPFM(T)
        PWM = learnPWM(PFM)
        col_i = get_column(PWM, i)
        # Determine the consensus based C_i
        C_i = col_i.index(max(col_i))
        # Calculate dependence S_i between C_i and other positions 


   
    return 0

if __name__ == '__main__':
    
    ### input data
    train_real_file = 'Data/hw3_train_real'
    train_false_file = 'Data/hw3_train_false'
    test_file = 'Data/hw3_test_real'
    #test_file = 'Data/hw3_test_false'
    sequences_real = read_sequences(train_real_file)
    sequences_false = read_sequences(train_false_file)
    sequences_test = read_sequences(test_file)
    PWM_positive = []
    PWM_negative = []
    ### Learn PWM_positive
    #PFM_positive = learnPFM(sequences_real)
    #PWM_positive = learnPWM(PFM_positive)
    #doctest.testmod()
    ### Learn PWM_negative
    #PFM_negative = learnPFM(sequences_false)
    #PWM_negative = learnPWM(PFM_negative)

    #####  MDD algorithm #####
    ### Init T, P
    T = sequences_real
    ncol = len(T[0][0])
    P = range(ncol)
    ### Built MDD Tree
    tree = find_MDD_subtree(T, P)
    
    eval_ChiSq([[4, 2], [0, 2], [0, 1], [0, 1]])


    doctest.testmod()
    


    print 'complied :D'
    

    
