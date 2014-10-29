#!/user/bin/env python
# Student: Duy Nguyen
# Solution to question ?,  homework ?
#
# You should be able to run the program as of the following commands


import sys
import io
import getopt

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

if __name__ == '__main__':
    read_score('score.txt')
#    print(read_score('score.txt'))

