#!/user/bin/env python
# Author: Duy Nguyen
# BMI/CS_776, Homework 3, Question 1
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
    """ Learn PFM (position frequency matrix) from positive examples
    
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

# def eval_prob_Zp(seq_x, pos_j, lengthW, lengthL, PWD):
#     """ Compute the P(X_i  | Z_ij = 1, p); pos_j is the real index, not python index (i.e. in python , start with 0)
#     >>> eval_prob_Zp('GCTGTAG', 3, 3, 7, [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]])
#     7.812500000000002e-06
#     >>> eval_prob_Zp('GCTGTAG', 1, 3, 7, [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]])
#     2.34375e-05
#     >>> eval_prob_Zp('GCTGTAG', 2, 3, 7, [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]])
#     0.00018750000000000003
#     """
#     index = {'A':0, 'C':1, 'G':2, 'T':3}
#     #seq_x = sequences[index_seq][0]
#     # Backgroud probability vector p0
#     p0 = get_column(PWD, 0)
#     # Compute 'before motif'
#     before = 1
#     if (pos_j > 1):
#         # Warning: k is 1 unit less than the 'real' position
#         for k in range(pos_j - 1):
#             p_c0 = p0[index[ seq_x[k] ] ]
#             before *= p_c0
#     if (pos_j < 1):
#         print "Starting position of the motif cannot be < 1"
#     # Compute 'motif'            
#     motif = 1
#     for k in range(pos_j - 1, pos_j + lengthW - 1):
#         col_k = get_column(PWD, k + 1 - pos_j + 1)
#         index_k = index[seq_x[k] ]
#         motif *= col_k[index_k]
#     # Compute 'after motif'            
#     #import pdb; pdb.set_trace()
#     after = 1
#     for k in range(pos_j + lengthW - 1, lengthL):
#         p_c0 = p0[index[ seq_x[k] ] ]
#         after *= p_c0
#     return before * motif * after    

def eval_prob_Zp(seq_x, pos_j, lengthW, lengthL, PWD):
    """ Compute the P(X_i  | Z_ij = 1, p); pos_j is the real index, not python index (i.e. in python , start with 0)
    """
    index = {'A':0, 'C':1, 'G':2, 'T':3}
    #seq_x = sequences[index_seq][0]
    if (pos_j > 1):
        # Warning: k is 1 unit less than the 'real' position
        for k in range(pos_j - 1):
            p_c0 = p0[index[ seq_x[k] ] ]
            before *= p_c0
    if (pos_j < 1):
        print "Starting position of the motif cannot be < 1"
    # Compute 'motif'            
    #import pdb; pdb.set_trace()
    motif = 1
    for k in range(pos_j - 1, pos_j + lengthW - 1):
        col_k = get_column(PWD, k)
        index_k = index[seq_x[k] ]
        motif *= col_k[index_k]
    # Compute 'after motif'            
    
    return  motif    

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

if __name__ == '__main__':
    
    ### input data
    train_real_file = 'Data/hw3_train_real'
    train_false_file = 'Data/hw3_train_false'
    #test_file = 'Data/hw3_test_real'
    test_file = 'Data/hw3_test_false'
    sequences_real = read_sequences(train_real_file)
    sequences_false = read_sequences(train_false_file)
    sequences_test = read_sequences(test_file)
    PWM_positive = []
    PWM_negative = []
    ### Learn PWM_positive
    PFM_positive = learnPFM(sequences_real)
    PWM_positive = learnPWM(PFM_positive)
    #doctest.testmod()
    ### Learn PWM_negative
    PFM_negative = learnPFM(sequences_false)
    PWM_negative = learnPWM(PFM_negative)
    ### Compute score(sequence) for sequence in test_file
    scores = []
    scores = evalScore(sequences_test, PWM_positive, PWM_negative )
    print scores[:10]
#PFM =  [[3, 6, 1, 0, 0, 6, 7, 2, 1], [2, 2, 1, 0, 0, 2, 1, 1, 2], [1, 1, 7, 10, 0, 1, 1, 5, 1], [4, 1, 1, 0, 10, 1, 1, 2, 6]]

    print 'complied :D'
   

    
