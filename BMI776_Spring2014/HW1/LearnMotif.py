#!/user/bin/env python
# Author: Duy Nguyen
# BMI/CS_776, Homework 1
# MEME OOPS Algorithm
# You should be able to run the program as of the following command
# python LearnMotif.py sequences file width model file positions file

import sys
import io
import getopt
import numpy as np


def read_sequences(filename):
    " Read sequences from file test01.txt "
    sequences = []
    with open(filename, 'r') as file:
        for line in file:
            seq = line.strip().split()
            sequences.append(seq)
    return sequences

def get_column(mat, i):
    " get the ith column of a given matrix"
    return [row[i] for row in mat]

def init_PWD(lengthW):
    " Initialize PWD matrix from the uniform [0,1] "
    mat = [[0 for x in range(lengthW+1)] for y in range(4)]    
    #import pdb; pdb.set_trace()
    for i in range(lengthW + 1):
        prob = np.random.randint(1, 100, 4)
        tot = sum(prob)
        for j in range(3):
           mat[j][i] = round(float(prob[j]) / tot, 3)
        col = [row[i] for row in mat]
        mat[3][i] = 1 - sum(col)   
    return mat

def eval_product(list):
    " Compute the product of a list"
    p = 1
    for i in list:
        p *= i
    return p

def eval_prob_Zp(seq_x, pos_j, lengthW, lengthL, PWD):
    " Compute the P(X_i  | Z_ij = 1, p)"
    index = {'A':0, 'C':1, 'G':2, 'T':3}
    #seq_x = sequences[index_seq][0]
    # Backgroud probability vector p0
    p0 = get_column(PWD, 0)
    # Compute 'before motif'
    before = 1
    if (pos_j > 1):
        # Warning: k is 1 unit less than the 'real' position
        for k in range(pos_j - 1):
            p_c0 = p0[index[ seq_x[k] ] ]
            before *= p_c0
    if (pos_j < 1):
        print "Starting position of the motif cannot be < 1"
    # Compute 'motif'            
    motif = 1
    for k in range(pos_j - 1, pos_j + lengthW - 1):
        col_k = get_column(PWD, k + 1 - pos_j + 1)
        index_k = index[seq_x[k] ]
        motif *= col_k[index_k]
    # Compute 'after motif'            
    #import pdb; pdb.set_trace()
    after = 1
    for k in range(pos_j + lengthW - 1, lengthL):
        p_c0 = p0[index[ seq_x[k] ] ]
        after *= p_c0
    return before * motif * after    


if __name__ == '__main__':
    sequences = []
    lengthW = 3
    lengthL = 7
    # create  the probability weight matrix filled with 0: rows' orders = a, c, g, t
    PWD = [[0 for x in range(lengthW + 1)] for y in range(4)]    
    # Initialize PWD
    #PWD = init_PWD(lengthW)
    PWD = [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]]

    sequences = read_sequences('Data/test01.txt')
    seq_x = sequences[0][0]
    foo = eval_prob_Zp('GCTGTAG', 3, lengthW, lengthL, PWD)
    print foo
    



    

