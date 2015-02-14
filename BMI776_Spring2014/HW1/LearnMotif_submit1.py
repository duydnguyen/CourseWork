#!/user/bin/env python
# Author: Duy Nguyen
# BMI/CS_776, Homework 1
# MEME OOPS Algorithm
# You should be able to run the program as of the following command
# python LearnMotif.py sequences file width model file positions file num_seed max_iteration epsilon
# Example: 
# python LearnMotif.py hw1_hidden_motif.txt 14 model_file positions_file 100 10000 0.001

import sys
import io
import getopt
import numpy as np
import doctest
from math import log

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

def init_PWD_seed(lengthW, seed):
    " Initialize PWD matrix from the random number with the given seed "
    mat = [[0 for x in range(lengthW+1)] for y in range(4)]    
    np.random.seed(seed)
    prob_main = np.random.randint(1, 10000, 4 * (lengthW + 1))
    #import pdb; pdb.set_trace()
    for i in range(lengthW + 1):
        start = 4 * i 
        end = start + 4
        prob = prob_main[start:end]
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
    """ Compute the P(X_i  | Z_ij = 1, p); pos_j is the real index, not python index (i.e. in python , start with 0)

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

def E_step_update(sequences, lengthN, lengthW, lengthL, PWD):
    """ Perform the E-step to estimate matrix Z(t). Z(t) must have row sum to 1.
    This function returns matrix matZ for Z(t) and also the log-likelihood logP(D|matP(t-1))
    """
    logL = 0
    tot_log = 0
    matZ = init_mat(lengthN, lengthL - lengthW + 1)
    for i in range(lengthN):
        seq_x = sequences[i][0]
        sum_P = 0
        for j in range(1, lengthL - lengthW + 2):
            #import pdb; pdb.set_trace()
            matZ[i][j-1]=  eval_prob_Zp(seq_x, j, lengthW, lengthL, PWD)
            sum_P += eval_prob_Zp(seq_x, j, lengthW, lengthL, PWD)
        matZ[i] = normalize(matZ[i])
        sum_P = log(sum_P)
        tot_log += sum_P
    logL = -lengthN * log(lengthL - lengthW + 1) + tot_log
    return matZ, logL    


def find_index(seq_x, char, index_k, lengthW, lengthL ):
    """ Given seq X, find the set of indices j <= L-W+1 such that X_{j+1-1} = char
    note that j indexed in MEME (start at 1, not 0)
    j in {1, 2, ..., L-W+1}; k in {1,2,..., W}

    >>> find_index('ACAGCA', 'A', 1, 3, 6)
    [1, 3]
    >>> find_index('AGGCAG', 'A', 1, 3, 6)
    [1]
    >>> find_index('TCAGTC', 'A', 1, 3, 6)
    [3]
    
    >>> find_index('ACAGCA', 'A', 2, 3, 6)
    [2]
    >>> find_index('AGGCAG', 'A', 2, 3, 6)
    [4]
    >>> find_index('TCAGTC', 'A', 2, 3, 6)
    [2]

    >>> find_index('ACAGCA', 'A', 3, 3, 6)
    [1, 4]
    >>> find_index('AGGCAG', 'A', 3, 3, 6)
    [3]
    >>> find_index('TCAGTC', 'A', 3, 3, 6)
    [1]

    >>> find_index('ACAGCA', 'C', 1, 3, 6)
    [2]
    >>> find_index('AGGCAG', 'C', 1, 3, 6)
    [4]
    >>> find_index('TCAGTC', 'C', 1, 3, 6)
    [2]

    >>> find_index('ACAGCA', 'C', 2, 3, 6)
    [1, 4]
    >>> find_index('AGGCAG', 'C', 2, 3, 6)
    [3]
    >>> find_index('TCAGTC', 'C', 2, 3, 6)
    [1]

    >>> find_index('ACAGCA', 'C', 3, 3, 6)
    [3]
    >>> find_index('AGGCAG', 'C', 3, 3, 6)
    [2]
    >>> find_index('TCAGTC', 'C', 3, 3, 6)
    [4]
    """
    index_j = []
    for j in range(1, lengthL - lengthW + 1 + 1):
        if seq_x[j + index_k - 1 - 1] == char:
            index_j.append(j)
    return index_j
                   

def total_char(seq_x, char):
    " Given seq X, count the total # of c's"
    total = 0
    for i in range(len(seq_x)):
        if seq_x[i] == char:
            total += 1
    return total
    

def eval_n_ck(sequences, lengthN, lengthW, lengthL, matZ):
    """ Evaluate matrix matN = [ n_ck ] where n_ck = number of the char 'c' occuring at position k over all seqs  

    >>> eval_n_ck([['ACAGCA'], ['AGGCAG'], ['TCAGTC']], 3, 3, 6, [[0.1, 0.7, 0.1, 0.1], [0.4, 0.1, 0.1, 0.4], [0.2, 0.6, 0.1, 0.1]])
    [[3.0999999999999996, 0.7000000000000001, 1.7000000000000002, 0.5], [2.5, 1.7000000000000002, 0.5, 0.30000000000000004], [1.7999999999999998, 0.4, 0.7, 2.1], [2, 0, 0, 0]]

    """
    # matN is a 4-by-(W + 1) matrix
    matN = init_mat(4, lengthW + 1)
    index = {'A':0, 'C':1, 'G':2, 'T':3}
    dic = {0:'A', 1:'C', 2:'G', 3:'T'}
    # index i: 0->A, 1->C, 2->G, 3->T
    # import pdb; pdb.set_trace()
    ## CASE k > 0
    for i in range(3):
        # k in {1, ..., W}; k = 0 will compute later
        for k in range(1, lengthW + 1):
            for seq_index in range(lengthN):
                seq_x = sequences[seq_index][0]
                char = dic[i]
                index_j = find_index(seq_x, char, k, lengthW, lengthL)
                for jj in index_j:
                    matN[i][k] += matZ[seq_index][jj-1]

    ## CASE K = 0
    # Compute n_c = total of char c in the data set
    n_c = []
    # index i for {A, C, G, T}
    for i in range(4):
        char = dic[i]
        total = 0
        for seq_index in range(lengthN):
            seq_x = sequences[seq_index][0]
            total += total_char(seq_x, char)
        n_c.append(total)
        matN[i][0] = n_c[i] - sum(matN[i])

    return matN 

def M_step(sequences, lengthN, lengthW, lengthL, matZ):
    """ Perform the M-step to estimate matrix PWD
    
    >>> M_step([['ACAGCA'], ['AGGCAG'], ['TCAGTC']], 3, 3, 6, [[0.1, 0.7, 0.1, 0.1], [0.4, 0.1, 0.1, 0.4], [0.2, 0.6, 0.1, 0.1]])
    [[0.30597014925373134, 0.25, 0.391304347826087, 0.21739130434782608], [0.2611940298507463, 0.39705882352941174, 0.21739130434782608, 0.18840579710144928], [0.20895522388059704, 0.20588235294117643, 0.24637681159420288, 0.4492753623188406], [0.2238805970149254, 0.14705882352941174, 0.14492753623188406, 0.14492753623188406]]

    """
    matP = init_mat(4, lengthW + 1)
    matN = eval_n_ck(sequences, lengthN, lengthW, lengthL, matZ)
    # index k for column
    for k in range(lengthW + 1):
        vec = get_column(matN, k)
        total = sum(vec)
        # index i for row
        for i in range(4):
            matP[i][k] = float(vec[i] + 1)/ (total + 4 )
    return matP

def eval_LogL(sequences, lengthN, lengthW, lengthL, matP):
    "Given the PWD matrix matP, evaluate the log-likelihood log P(D | matP(t))"
    logL = 0
    tot_log = 0
    for i in range(lengthN):
        #import pdb; pdb.set_trace()    
        seq_x = sequences[i][0]
        sum_P = 0 
        for j in range(1, lengthL - lengthW + 2):
            sum_P += eval_prob_Zp(seq_x, j, lengthW, lengthL, matP)
        sum_P = log(sum_P)
        tot_log += sum_P    
    #import pdb; pdb.set_trace()    
    logL = -lengthN * log(lengthL - lengthW + 1) + tot_log
    return logL

def output_model(lengthW, matP, model_file):
    " Write the PWD matrix matP to file "
    f = open(model_file, 'w')
    for row in range(4):
        for col in range(lengthW + 1):
            f.write(str(matP[row][col])+'\t')
        f.write('\n')
    f.close()

def find_site(seq_index, matZ_best):
    " Output best motif starting position given the sequence (real index, i.e. starts at 1)"
    val = max(matZ_best[seq_index-1])
    for i in range(len(matZ_best[0])):
        if ( matZ_best[seq_index-1][i] < val + 0.1 ) and ( matZ_best[seq_index-1][i] > val - 0.1 ):
            return i

def output_positions(lengthN, matZ, positions_file):
    " Write positions_file which contains the best motif starting positions for each sequence "
    f = open(positions_file, 'w')
    for row in range(4):
        for col in range(lengthW + 1):
            f.write(str(matP[row][col])+'\t')
        f.write('\n')
    
    for i in range(1, lengthN + 1):
        pos = find_site(i, matZ)
        f.write('Sequence ' + str(i) + ' starts at ' + str(pos) + '\n')
    f.close()
    
def main(argv):
    """Take arguments from the following command:
    python LearnMotif.py sequences file width model file positions file num_seed max_iteration epsilon
    """
    try:
        opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # ## Error handlings: LATER!
    # if len(args) < 5:
    #     sys.exit('There should be a total of 5 arguments, at least one missing')
    # if len(args) > 5:
    #     sys.exit('There should be a total of 5 arguments, there is at least one extra argument')
    sequences_file = args[0]
    width = int(args[1])
    model_file = args[2]
    positions_file = args[3]
    num_seed = args[4]
    max_iteration = args[5]
    epsilon = args[6]
    return [sequences_file, width, model_file, positions_file, num_seed, max_iteration, epsilon]
    
if __name__ == '__main__':
    
    ### Inputs
    inputs = main(sys.argv[1:])
    sequences_file = inputs[0]
    lengthW = inputs[1]
    model_file = inputs[2]
    positions_file = inputs[3]
    tot_seed = inputs[4]
    max_t = inputs[5]
    sequences = []
    sequences = read_sequences(sequences_file)
    # Cut-off threshold for log-likelihood
    epsilon = inputs[6]
    ### Evaluate length of a sequence (lengthL) and number of sequences (lengthN)
    lengthL = len(sequences[0][0])
    lengthN = len(sequences)


    #model_file = 'model_file'
    #positions_file = 'positions_file'
    #sequences = []
    #sequences = read_sequences('Data/hw1_hidden_motif.txt')
    #lengthW = 14
    #lengthL = 200
    #lengthN = len(sequences)
    ### seed =  number of starting points
    tot_seed = 5
    count_seed = 0
    logL_seed = float('-inf')
    best_seed = -1
    matZ_best = []
    matP_best = []
    while (count_seed < tot_seed):
        count_seed += 1
        print '\n'
        print '*****************************'
        print '+++ Starting point number = ' + str(count_seed)
        ### Initilize PWD matrix
        s = 2 * (count_seed - 1) + 1 
        matP = init_PWD_seed(lengthW, seed = s )
        ### Run first iteration: logL_prev from matP(t-1)
        matZ, logL_prev = E_step_update(sequences, lengthN, lengthW, lengthL, matP)
        matP = M_step(sequences, lengthN, lengthW, lengthL, matZ)
        ### index for interation
        t = 0
        max_t = 10000
        ### Cut-off threshold for log-likelihood
        epsilon = 0.001
        check = False # False if change in logL >= epsilon
        while ( (t <= max_t) and (check == False) ):
            t += 1
            # E-step: re-estimate Z(t), and compute logL from the new matP(t)
            matZ, logL = E_step_update(sequences, lengthN, lengthW, lengthL, matP)
            # M-step: re-estimate matP(t+1)
            matP = M_step(sequences, lengthN, lengthW, lengthL, matZ)
            # Check for change in logL < epsilon
            if abs(logL - logL_prev) < epsilon:
                check = True
                print logL
                print '++++++++ Current best seed = ' + str(best_seed) + ' with logL = ' + str(logL_seed)
                # This stores the best logL_seed = best logL with all seeds considered so far
                if logL > logL_seed:
                    logL_seed = logL
                    best_seed = s
                    matZ_best = matZ
                    matP_best = matP
                    print '++++++++ Current best seed = ' + str(best_seed) 
            else:
                logL_prev = logL
    ## Print the optimal logL
    print '+++ Optimal logL' + str(logL_seed)
    print '+++ Best seed = ' + str(best_seed)

    ## Output: model_file which contains PWD matrix
    output_model(lengthW, matP_best, model_file)
    ## Output: positions_file which contains the best motif starting positions for each sequence
    output_positions(lengthN, matZ_best, positions_file)
