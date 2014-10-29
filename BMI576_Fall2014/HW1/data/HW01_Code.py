#!/user/bin/env python
# Student: Duy Nguyen
# Solution to question 1,  homework 1
#
# Read FASTA file
# You should be able to run the program as of the folloing commands
# python HW01_Code.py yeast_chrom.fasta chr01 20 30 +

import sys
import io
import getopt

def read_FASTA_sequence(file):
    "Read the sequence given the chromosome"
    seq = ''
    for line in file:
        if not line or line[0] == '>':
            return seq
        seq += line[:-1]

    return seq

def get_id(line):
    "Get the chromosome number"
    return line.replace('>', '').replace('\n','')


def FASTA_search_by_id(id, file):
    for line in file:
        if (line[0] == '>' and str(id) == get_id(line)):
            return read_FASTA_sequence(file)


def search_FASTA_file_by_id(id, filename):
    "Search FASTA file and read only the seq with the corresponding chromosome"
    id = str(id)
    with open(filename) as file:
        return FASTA_search_by_id(id, file)

def eval_Seq(seq, start, end, strand = '+'):
    res = ''
    if strand == '+':
        res = seq[(start-1):end]
    else:
        res = rcomplement(seq[(start-1):end])
        
    return res

def rcomplement(s):
    "Find the reverse complement of the sequence"
    basecomplement = {'A': 'T', 'C':'G', 'G':'C', 'T':'A'}
    letters = list(s)
    # find the complement
    letters = [basecomplement[base] for base in letters]
    
    return ''.join(letters)[::-1]

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hg:d", ["help", "grammar="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    FileName = args[0]
    Chr = args[1]
    Start = int(args[2])
    End = int(args[3])
    Strand = args[4]

    ## Error Handlings                                                                                                 
    if Chr[:3] != "chr":
        sys.exit('Chromosome argument is not in right format')
    if Start < 1:
        sys.exit('Start Position should be greater than 0')
    if not (Strand == '+' or Strand == '-'):
        sys.exit('Input of Strand must be +/-')

    id = Chr
    seq = search_FASTA_file_by_id(id, FileName)
    if End > len(seq):
        sys.exit('End Position exceeds chromosome length')
    seq_frag = eval_Seq(seq, start = Start, end = End, strand = Strand)
    print(seq_frag)


if __name__ == '__main__':
    main(sys.argv[1:])
    
