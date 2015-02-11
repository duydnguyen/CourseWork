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
import doctest

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
    """ Compute the P(X_i  | Z_ij = 1, p); pos_j is the real index, not python index (i.e. start with 0)

    >>> eval_prob_Zp('GCTGTAG', 3, 3, 7, [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]])
    7.812500000000002e-06
    >>> eval_prob_Zp('GCTGTAG', 1, 3, 7, [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]])
    2.34375e-05
    >>> eval_prob_Zp('GCTGTAG', 2, 3, 7, [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]])
    0.00018750000000000003

    """
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

def init_mat(nrow, ncol):
    " Initialize nrow-by-ncol matrix with 0"
    return  [[0 for x in range(ncol)] for y in range(nrow)]  
    
def normalize(list):
    " Normalize a list"
    tot = sum(list)
    for i in range(len(list)):
        list[i] = float(list[i]) / tot
    return list

def E_step(sequences, lengthN, lengthW, lengthL, PWD):
    """ Perform the E-step to estimate matrix Z(t). Z(t) must have row sum to 1.
    
    """
    matZ = init_mat(lengthN, lengthL - lengthW + 1)
    for i in range(lengthN):
        seq_x = sequences[i][0]
        for j in range(1, lengthL - lengthW + 2):
            #import pdb; pdb.set_trace()
            matZ[i][j-1]=  eval_prob_Zp(seq_x, j, lengthW, lengthL, PWD)
        matZ[i] = normalize(matZ[i])
    return matZ    

def find_index(seq_x, char, index_k, lengthW, lengthL ):
    """ Given seq X, find the set of indices j <= L-W+1 such that X_{j+1-1} = char
    note that j indexed in python (start at 0, not 1)
    >>> find_index('ACAGCA', 'A', 1, 3, 6)
    [0, 2]
    >>> find_index('ACAGCA', 'C', 1, 3, 6)
    [1]
    >>> find_index('ACAGCA', 'G', 1, 3, 6)
    [3]
    >>> find_index('ACAGCA', 'T', 1, 3, 6)
    []
    >>> find_index('AGGCAG', 'A', 1, 3, 6)
    [0]
    >>> find_index('ACAGCA', 'C', 1, 2, 6)
    [1, 4]
    >>> find_index('ACAGCA', 'C', 3, 3, 6)
    [2]
    
    """
    #>>> find_index('ACAGCA', 'C', 4, 3, 6)
    seq = seq_x[0:(lengthL - lengthW + 1)]
    index_j = []
    for j in range(len(seq)):
        if seq_x[j + index_k - 1] == char:
            index_j.append(j)
    return index_j
                   

def eval_n_ck(sequences, lengthN, lengthW, lengthL, matZ):
    """ Evaluate matrix matN = [ n_ck ] where n_ck = number of the char 'c' occuring at position k over all seqs  

    """
    matN = init_mat(lengthW + 1, 4)
    
    k = 1
    for i in range(lengthN):
        seq_x = sequences[i][0]
        for j in range(1, lengthL - lengthW + 2):
            k = 1
            

    return matN 

def M_step():
    """ Perform the M-step to estimate matrix PWD
    """
    return 0

if __name__ == '__main__':
    sequences = []
    lengthW = 3
    lengthL = 7
    ## Create  the probability weight matrix filled with 0: rows' orders = a, c, g, t
    PWD = init_mat(lengthW + 1, 4)
    ## Initialize PWD
    #PWD = init_PWD(lengthW)
    PWD = [[0.25, 0.1, 0.5, 0.2], [0.25, 0.4, 0.2, 0.1], [0.25, 0.3, 0.1, 0.6], [0.25, 0.2, 0.2, 0.1]]
    #sequences = read_sequences('Data/test01.txt')
    ## number of DNA sequences
    lengthN = len(sequences)
    ## Create  the Z matrix of latent variables for motif positions filled with 0: 
    matZ = init_mat(lengthN, lengthL - lengthW + 1)
    #matZ = E_step(sequences, lengthN, lengthW, lengthL, PWD)
    #print matZ
    #seq_x = sequences[0][0]
    #foo = eval_prob_Zp('GCTGTAG', 5, lengthW, lengthL, PWD)
    #print foo
    #doctest.testmod()

    ### THIS IS FOR eval_n_ck()
    sequences = [['ACAGCA'], ['AGGCAG'], ['TCAGTC']]
    matZ = [[0.1, 0.7, 0.1, 0.1], [0.4, 0.1, 0.1, 0.4], [0.2, 0.6, 0.1, 0.1]]
    matN = eval_n_ck(sequences, 3, 3, 6, matZ)
    doctest.testmod()
    print matN












    
    



    

