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
            #import pdb; pdb.set_trace()
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
    foo = list(set(relations.keys() + relations.values()))
    return len(foo)

def inverseDict(Dict):
    "Inverse a dictionary, also for case of non-unique map. This is used for computing childrens of a given internal node"
    invDict = {}
    for k, v in Dict.iteritems():
        # This creates a list to store nodes of inverse map
        invDict[v] = invDict.get(v, [])
        invDict[v].append(k)
    return invDict

def addList(L1, L2):
    return [x + y for x, y in zip(L1, L2)]

def mapEval(relations):
    "Create the mapping from node_unique (map key) to {1,..., nodeNum} (map value)"
    node_unique = list(set(relations.keys() + relations.values()))        
    mapping = {}
    index = 1
    for x in node_unique:
        mapping[x] = index
        index += 1 
    return mapping

def relabel(relations, mapping):
    "Relabel nodes in var relations to avoid the non-consecutive labeling case."
    relations_new = {}
    for k, v in relations.iteritems():
        relations_new[mapping[k]] = mapping[v]
    return relations_new

def relabel_assign(assign, mapping):
    "Relabel leaf nodes in var assign"
    assign_new = {}
    for k, v in assign.iteritems():
        assign_new.update({mapping[k]:v})
    return assign_new

def costEval(score, relations, assign):
    "Evaluate the cost matrix by the weighted parsimony algorithm"
    import pdb; pdb.set_trace()
    ## Initialize, relabeling, and find root of tree
    leaves = []
    relations_new = {}
    assign_new = {}
    mapping = {}
    # Evaluate number of nodes
    nodeNum = nodeNum_Eval(relations)
    # Create the mapping for relabeling
    mapping = mapEval(relations)
    # Relabel nodes in relations
    relations_new = relabel(relations, mapping)
    # Relabel leaves in assign
    assign_new = relabel_assign(assign, mapping)
    # root of tree
    root = list( set(range(1,nodeNum+1)) - set(relations_new.keys()))[0]

    ## Initialize the Cost Matrix. Ordering as score matrix:  Col1 = 'a', Col2 = 'c', Col3 = 'g', Col4 = 't'
    index = {'a':1, 'c':2, 'g':3, 't':4}
    Cost = [[0 for x in range(4)] for y in range(nodeNum)]
    # Initialize the Cost(leaf)
    for x in assign_new:
        leaves.append(x)
        Cost[x-1] = [float("inf"), float("inf"), float("inf"), float("inf")]
        Cost[x-1][index[assign_new[x]]-1] = 0
    # Evaluate Cost matrix for internal nodes
    intNodes = [x for x in range(1,nodeNum+1) if not x in leaves]
    inv_relations = inverseDict(relations_new)
    for i in intNodes:
        child1, child2 = inv_relations[i][0], inv_relations[i][1]
        for j in range(4):
            # import pdb; pdb.set_trace()
            min1 = min( addList(Cost[child1 - 1], score[j]))
            min2 = min( addList(Cost[child2 - 1], score[j]))
            # # for debugging
            # print "+++ parent node=%s, child1=%s, child2=%s" % (i, child1, child2)
            # print("+++ score[j] = ", score[j])
            # print("+++ Cost[child1 - 1] = ", Cost[child1-1])
            # print("+++ Cost[child2 - 1] = ", Cost[child2-1])
            # #
            Cost[i-1][j] = min1 + min2
    print(Cost)
    print "Cost of the tree = %s" % min(Cost[nodeNum-1])
    return Cost

    

if __name__ == '__main__':
    score = read_score('Tests/Q1_Test5/score.txt')
    relations = read_tree('Tests/Q1_Test5/tree.txt')
    print(relations)
    assign = read_assign('Tests/Q1_Test5/assign.txt')
    print(assign)
    costEval(score, relations, assign)



