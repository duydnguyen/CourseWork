#!/user/bin/env python
# Author: Duy Nguyen
# BMI/CS_776, 
# Maximum Dependence Decomposition (MDD) Model
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

def eval_Matches_consensus(col_i, col_j, char, C_i):
    """ Given a char in {A,C,G,T} and consensus base C_i, evaluate matches/non-matches between two sequence columns
    >>> eval_Matches_consensus('GTTTGCTA', 'TTAAAATA', 'A', 'T')
    [2, 3]

    >>> eval_Matches_consensus('GTTTGCTA', 'TTAAAATA', 'T', 'T')
    [2, 1]

    >>> eval_Matches_consensus('GTTTGCTA', 'TTGTTGCC', 'T', 'T')
    [2, 2]

    >>> eval_Matches_consensus('GTTTGCTA', 'TTGTTGCC', 'C', 'T')
    [1, 1]
    """
    pos = find_Indices(col_j, char)
    combined = []
    for i in range(len(pos)):
        s = col_i[pos[i]] + col_j[pos[i]]
        combined.append(s)
    # eval the number of matches
    matches = 0
    #check = [False for i in range(len(combined))]
    pos_Ci = find_Indices(col_i, C_i)
    pos_common = set(pos) & set(pos_Ci)
    matches = len(pos_common)
    non_matches = len(pos) - matches
    return [matches, non_matches]


def eval_ChiSq(table):
    """ Given a 4-by-2 contigency table, evaluate the Chi Square Statistics
    >>> eval_ChiSq([[4, 2], [0, 2], [0, 1], [0, 1]])
    4.444444444444445

    >>> eval_ChiSq([[1, 2], [3, 4], [5, 6], [7, 8]])
    0.19168831168831166
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
            if (table_expected[i][j] != 0):
                Chi += float(table[i][j] - table_expected[i][j] ) ** 2 / table_expected[i][j]
    return Chi

def eval_Table(sequences, index_i, index_j, C_i):
    """ Given set of sequences and indices i, j , eval a 4-by-2 contigency table
    
    >>> eval_Table( [['GAGGTAAAC'],['TCCGTAAGT'],['CAGGTTGGA'],['ACAGTCAGT'],['TAGGTCATT'],['TAGGTACTG'],['ATGGTAACT'],['CAGGTATAC'],['TGTGTGAGT'],['AAGGTAAGT']], 0, 1, 'T')
    [[2, 4], [1, 1], [1, 0], [0, 1]]
    
    """
    
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    table =  [[0 for x in range(2)] for y in range(4)]    
    col_i = extract_col(sequences, index_i )
    col_j = extract_col(sequences, index_j )

    for i in range(4):
        char = index[i]
        table[i] = eval_Matches_consensus(col_i, col_j, char, C_i)
    return table

def eval_Si(sequences, index_i, C_i):
    " Given set of sequences, index i and consensus base C_i, compute S_i = \sum_{j\neq i} ChiSq(C_i, x_j) "
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    ncol = len(sequences[0][0])
    table =  [[0 for x in range(2)] for y in range(4)]    
    col_i = extract_col(sequences, index_i )
    Si = 0
    for j in range(ncol):
        if (j != index_i):
            table = eval_Table(sequences, index_i, j, C_i)
            Si += eval_ChiSq(table) 
    return Si
        
def eval_Si_store(T, P, PWM):
    ' Eval Si_store[i] for all i in P'
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    Si_store = []
    Si_store_full = []
    ncol = len(T[0][0])
    for i in range(ncol):
        #PFM = learnPFM(T)
        col_i = get_column(PWM, i)
        # Determine the consensus base C_i: get the nucleotide with max probability in col_i
        C_i = index[col_i.index(max(col_i))]
        #consensus.append(C_i)
        # Calculate dependence S_i between C_i and other positions 
        Si_store.append(eval_Si(T, i, C_i))
        Si_store_full.append(eval_Si(T, i, C_i))
    
    #print '~~~~~~ NEW eval_Si_store() ~~~~'
    #print '+++ P = % s' % P
    #print '+++ Si_store = % s' % Si_store 
    diff =  [x for x in range(ncol) if x not in set(P)]
    #print '+++ diff = % s' % diff
    Si_store =  [item for i,item in enumerate(Si_store) if i not in diff]
    maxS = max(Si_store)
    # find i_max in full P i.e., range(9)
    i_max = Si_store_full.index(maxS)
    
    P_update = removeList(P, i_max)
    return Si_store_full, maxS, i_max, P_update

def split_Sequences(T, i_max, C_i):
    "Split sequences in to D_i+ and D_i- by consensus base C_i in MMD algorithm"
    Di_plus = []
    Di_minus = []
    for k in range(len(T)):
        seq = T[k][0]
        if (seq[i_max] == C_i):
            Di_plus.append([seq])
        else:
            Di_minus.append([seq])
    
    return Di_plus, Di_minus

def removeList(L, i):
    'Given i \in L, remove the i element from L'
    L_removed = [x for x in L if x != i]
    return L_removed
    
def find_MDD_subtree(T, P):
    " Find the tree using the Maximal Dependence Decomposition (MDD) algorithm"
    print '\n \n ~~~~~~~~~~~~~~~~~~~~~~~~~~~NEW RECURSION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    global Tree
    global Nodes
    global Store
    global Parent
    global Parent_Prob
    global Leaf_Prob
    global Leaf_Prob_index    
    global Check_Subtree
    ncol = len(T[0][0])
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    negate_index = {'G':'H', 'A':'B', 'T':'V', 'C':'D'}
    cutoff_seq = 399
    cutoff_ChiSq = 16.3
    ## Compute PWM given current T and P
    PFM = learnPFM(T)
    PWM = learnPWM(PFM)
    print '+++ current P = % s' % P
    #print '+++ current len(T) = % s' % len(T)
    #print '+++ current Nodes = % s' % Nodes

    ### STEP 1
    Si_store = []
    Si_store_foo, maxS, i_max_foo, P_update_foo = eval_Si_store(T, P, PWM)
    print '++++ the new i_max_foo  : % s and maxS = % s' % (i_max_foo, maxS)
    print '++++ the P_update_foo after removing i_max : % s' % P_update_foo
    
    ### STEP 2
    if (len(T) > cutoff_seq) and (maxS > cutoff_ChiSq):
        ### STEP 1
        Si_store = []
        Si_store, maxS, i_max, P_update = eval_Si_store(T, P, PWM)
        print '++++ the new i_max  : % s and maxS = % s' % (i_max, maxS)
        print '+++++ the P_update after removing i_max : % s' % P_update

        ## Choose the value i such that S_i is maximal
        #i_max = Si_store.index(maxS)
        
        ## Make a node with consensus base C_i as the test
        # and create a single-column PWM col_i for position i
        col_i = get_column(PWM, i_max)
        C_i = index[col_i.index(max(col_i))]
        ## Partition set of sequences T in to D_i+ and D_i-
        Di_plus, Di_minus = split_Sequences(T, i_max, C_i)
        print '+++ len(Di_plus) = % s ' % len(Di_plus)
        print '+++ len(Di_minus) = % s ' % len(Di_minus)
        ## Build Tree
        #Tree.append([Nodes, Nodes+1, Nodes +2])
        Store[Nodes + 1] = [i_max, C_i]
        Store[Nodes + 2] = [i_max, negate_index[C_i] ]
        Tree.append([Parent, Nodes+1, Nodes +2])
        Parent_Prob[Parent] = col_i
        Nodes += 2
        #print '+++ Nodes = % s' % Nodes
        print '+++ Store = % s' % Store
        ## left subtree
        #import pdb; pdb.set_trace()
        #diff = [i_max]
        # remove i_max element for P
        
        # print '++++ \n \n FOUND ERROR'
        # print 'i_max = % s' % i_max
        # print '+++ current P = % s' % P
        #P.remove(i_max)
        
        print '++++++++ RUNNING LEFT SUBTREE'
        Parent = Nodes - 1
        Check_Subtree = True
        print '+++ This is the P for left subtree P = % s' % P_update
        print '+++ This is the i_max left subtree i_max = % s' % i_max
        find_MDD_subtree(Di_plus, P_update)
        ## right subtree
        Parent = Nodes 
        print '++++++++ RUNNING RIGHT SUBTREE'
        Check_Subtree = False
        print '+++ \t This is the P for right subtree = % s' % P_update
        print '+++ This is the i_max right subtree i_max = % s' % i_max
        find_MDD_subtree(Di_minus, P_update)
    else:
        print '\n \n @@@@@@@@@@ THIS IS WHEN STOPPING CRIT. MET!!!'
        print '+++ current P = % s' % P
        #print 'Current Parent = % s' % Parent
        print '+++ current len(T) = % s' % len(T)
        print '++ Check_Subtree = % s' % Check_Subtree
        print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        
        diff =  [x for x in range(ncol) if x not in set(P)]
        print 'len of diff = % s' % diff 
        if (Check_Subtree):
            node_Label = len(diff) * 2
            Leaf_Prob[node_Label] = PWM
            Leaf_Prob_index[node_Label] = P
        else:
            node_Label = len(diff) * 2 + 1
            Leaf_Prob[node_Label] = PWM
            Leaf_Prob_index[node_Label] = P
        #print '+++ current Leaf_Prob = % s' % Leaf_Prob
    return 0

def built_MDDmodel(T):
    "Given a set of sequences (positive or negative), outout MDD model"
    global Tree
    global Nodes
    global Store
    global Parent
    global Parent_Prob
    global Leaf_Prob
    global Leaf_Prob_index
    global Check_Subtree
    global num_Nodes
    global ncol
    global Tree_struct
    ### Initilize all variables
    ncol = len(T[0][0])
    #P = range(ncol)
    Tree = []
    Nodes = 1 
    Store = {}
    Tree_struct = {}
    Parent = 1
    num_Nodes = 0
    # store a column of probabilty for parent (internal) nodes
    Parent_Prob = {}
    Leaf_Prob = {}
    Leaf_Prob_index = {}
    # check current states whether in left subtree or right subtree: Check_Subtree = True of on right subtree; False otherwise
    Check_Subtree = True
    find_MDD_subtree(T, range(ncol))
    num_Nodes = len(Store) + 1
    # check if a node is internal node or leaf: True = internal node; False = leaf
    Tree_struct = [False for i in range( num_Nodes+1)]
    print '\n \n +++++++ RESULTS +++++++'
    print '\n Store = % s , \n Tree = % s, \n Number of Nodes = % s' % (Store, Tree, num_Nodes)
    print '\n Parent Probabilty  = % s' % Parent_Prob
    #print '\n Leaf Probability  = % s' % Leaf_Prob
    print '\n Leaf Probability index = % s' % Leaf_Prob_index

    ### Find internal nodes of MDD Tree
    for i in range(len(Tree)):
        int_node = Tree[i][0]
        if len(Tree[i]) > 1:
            Tree_struct[int_node] = True
    #print '\n Tree_struct = % s' % Tree_struct
    print '\n Labels of nodes in MDD Tree  = % s' % Tree_struct
    return 0

def inverseDict(Dict):
    "Inverse a dictionary, also for case of non-unique map. This is used for computing childrens of a given internal node"
    invDict = {}
    for k, v in Dict.iteritems():
        # This creates a list to store nodes of inverse map
        invDict[v] = invDict.get(v, [])
        invDict[v].append(k)
    return invDict

def findParent(Tree):
    """ Given Tree, create a dictionary to store parent for each node
    >>> findParent([[1, 2, 3], [2, 4, 5], [4, 6, 7]])
    {2: 1, 3: 1, 4: 2, 5: 2, 6: 4, 7: 4}
    
    >>> findParent([[1, 2, 3], [3, 4, 5], [5, 6, 7]])
    {2: 1, 3: 1, 4: 3, 5: 3, 6: 5, 7: 5}
    """
    dict_Parent = {}
    for L in Tree:
        par = L[0]
        dict_Parent[ L[1] ] = par
        dict_Parent[ L[2] ] = par
    return dict_Parent

def findTreeDict(Tree):
    "Given Tree, convert it to dictionary"
    TreeDict ={}
    for L in Tree:
        par = L[0]
        TreeDict[par] = [L[1], L[2]]
    return TreeDict

def evalProb_leaf(sequence, leaf_index, leaf_PWM):
    "Evaluate probility of leaf node"
    prob = 1
    inverse_index = {'A':0, 'C':1, 'G':2, 'T':3}
    #import pdb; pdb.set_trace()
    for i in leaf_index:
        char_index = inverse_index[ sequence[i]]
        col = get_column(leaf_PWM, i)
        prob *= col[ char_index ] 
    return prob

def output_Prob(sequence, Node):
    'Eval the probability of a sequence given the current MDD model using RECURSION'
    global prob_seq 
    index = {0:'A', 1:'C', 2:'G', 3:'T'}
    inverse_index = {'A':0, 'C':1, 'G':2, 'T':3}
    negate_index = {'G':'H', 'A':'B', 'T':'V', 'C':'D'}
    dict_Parent = {}
    TreeDict = {}
    ## Find parent of nodes given Tree
    dict_Parent = findParent(Tree)
    TreeDict = findTreeDict(Tree)
    ## Case: Node is leaf
    if Tree_struct[Node] == False:
        leaf_index = Leaf_Prob_index[Node]
        leaf_PWM = Leaf_Prob[Node]
        print '+++ node is LEAF: leaf_index = % s' % leaf_index
        prob_leaf = evalProb_leaf(sequence, leaf_index, leaf_PWM)
        prob_seq *= prob_leaf
    ## Case: Node is internal
    else:
        child_left = TreeDict[Node][0]
        child_right = TreeDict[Node][1]
        i_max = Store[child_left][0] # 7
        Ci = Store[child_left][1] # G
        # x_pos =  base of sequence at i_max
        x_pos = sequence[i_max] 
        col_imax = Parent_Prob[Node]
        prob_seq *= col_imax[inverse_index[ x_pos ]]
        #import pdb; pdb.set_trace()
        print 'Node = % s, i_max = % s, Ci = % s, x_pos = % s, \n prob_seq = % s' % (Node, i_max, Ci, x_pos, prob_seq)
        if x_pos != Ci:
            print '\n @@@ Running child_right = % s ' % child_right
            output_Prob(sequence, child_right)
            print '+++ current prob_seq = % s' % prob_seq
        else:
            print '\n @@@ Running child_left = % s ' % child_left
            output_Prob(sequence, child_left)
            print '+++ current prob_seq = % s' % prob_seq
     
    return 0

def eval_score(MDD_p, MDD_n):
    'Eval score(sequence) = log ( P(seq|MDD+) / P(seq|MDD-)  )'
    scores = []
    for i in range(len(MDD_p)):
        scores.append( log(MDD_p[i]) - log(MDD_n[i])   )
        
    return scores
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
    
    ### Built MDD Tree
    ## 1. Built MDD_p Tree
    T = sequences_real
    built_MDDmodel(T)
    ##############################
    ## Compute prob(seq| MDD_p)
    prob_PositiveMDD = []
    for seq in sequences_test:
        seq_test = seq[0]
        prob_seq = 1        
        output_Prob(seq_test, 1)
        #print 'seq = % s, prob_seq = % s' % (seq_test, prob_seq)
        prob_PositiveMDD.append(prob_seq)

    ## 2. Built MDD_n Tree
    T = sequences_false
    built_MDDmodel(T)
    ##############################
    ## Compute prob(seq| MDD_n)
    prob_NegativeMDD = []
    for seq in sequences_test:
        seq_test = seq[0]
        prob_seq = 1        
        output_Prob(seq_test, 1)
        #print 'seq = % s, prob_seq = % s' % (seq_test, prob_seq)
        prob_NegativeMDD.append(prob_seq)

    #doctest.testmod()
    scores = eval_score(prob_PositiveMDD, prob_NegativeMDD)
    print scores
    print 'complied :D'
    

    
