#!/user/bin/env python
# Student: Duy Nguyen
# Solution to question 1,  homework 3
#
# You should be able to run the program as of the following commands


import sys
import io
import getopt
import operator

def read_score(filename):
    "Read score matrix from score.txt"
    # Initialize 2D list with dim = 4*4
    score = [[0 for x in range(4)] for y in range(4)]
    row = 0
    with open(filename, 'r') as file:
        file.readline() # omit the first line
        for line in file:
            vals = line.strip().split()
            score[row] = [float(vals[1]), float(vals[2]), float(vals[3]), float(vals[4])]
            row += 1
    return score

def read_tree(filename):
    "Read child/parent from tree.txt"
    # Initialize the dictionary to store child/parent relations
    relations = {}
    child = 0
    parent = 0
    with open(filename, 'r') as file:
        for line in file:
            vals = line.strip().split()
            child = int(vals[0].strip('S'))
            parent = int(vals[1].strip('S'))
            relations[child] = parent
    return relations

def read_assign(filename):
    "Read assignments at leaf nodes from assign.txt"
    # Initialize the dictionary for assignments: NO space format s1=a
    assign = {}
    leaf = 0
    character = ''
    with open(filename, 'r') as file:
        for line in file:
            vals = line.strip().split()
    for i in range(len(vals)):
        leaf = int(vals[i][1])
        character = vals[i][3]
        assign[leaf] = character
    return assign

def nodeNum_Eval(relations):
    "Evaluate number of nodes for the given tree"
    key,value = max(relations.iteritems(), key=lambda x:x[1])
    return value

def inverseDict(Dict):
    "Inverse a dictionary, also for case of non-unique map. This is used for computing childrens of a given internal node"
    invDict = {}
    for k, v in Dict.iteritems():
        # This creates a list to store nodes of inverse map
        invDict[v] = invDict.get(v, [])
        invDict[v].append(k)
    return invDict

def costEval(score, relations, assign):
    "Evaluate the cost matrix by the weighted parsimony algorithm"
    leaves = []
    # Evaluate number of nodes
    nodeNum = nodeNum_Eval(relations)
    # Initialize the Cost Matrix where Col1 = 'a', Col2 = 't', Col3 = 'g', Col4 = 'c'
    index = {'a':1, 't':2, 'g':3, 'c':4}
    Cost = [[0 for x in range(4)] for y in range(nodeNum)]
    # Initialize the Cost(leaf)
    for x in assign:
        leaves.append(x)
        Cost[x-1] = [float("inf"), float("inf"), float("inf"), float("inf")]
        Cost[x-1][index[assign[x]]-1] = 0
    # Evaluate Cost matrix for internal nodes
    intNodes = [x for x in range(1,nodeNum+1) if not x in leaves]
    inv_relations = inverseDict(relations)
    for i in intNodes:
        for j in range(3):
            child1, child2 = inv_relations[i][0], inv_relations[i][1]
            import pdb; pdb.set_trace()
            min1 = min( Cost[child1 - 1] + score[j])
            min2 = min( Cost[child2 - 1] + score[j])
            Cost[i-1][j] = min1 + min2
            #Cost[i-1][j] = min( Cost[child1 - 1] + score[j]) + min( Cost[child2 - 1] + score[j])  
    print(Cost)
    return Cost

    

if __name__ == '__main__':
    score = read_score('score.txt')
    relations = read_tree('tree.txt')
    assign = read_assign('assign.txt')
#    print(read_assign('assign.txt'))
#    print(read_tree('tree.txt'))
    
    costEval(score, relations, assign)



