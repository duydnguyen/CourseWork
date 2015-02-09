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

if __name__ == '__main__':
    sequences = []
    lengthW = 3
    # create  the probability weight matrix filled with 0: rows' orders = a, c, g, t
    PWD = [[0 for x in range(lengthW + 1)] for y in range(4)]    
    # Initialize PWD
    PWD = init_PWD(lengthW)

    sequences = read_sequences('Data/test01.txt')
    
    



    

