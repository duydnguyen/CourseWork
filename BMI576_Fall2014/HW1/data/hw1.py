#!/user/bin/env python
#
# Solution to question 1,  homework 1
#
# Read FASTA file

def read_FASTA_sequence(file):
    seq = ''
    for line in file:
        if not line or line[0] == '>':
            return seq
        seq += line[:-1]

    return seq

def get_id(line):
    return line.replace('>', '').replace('\n','')


def FASTA_search_by_id(id, file):
    for line in file:
        if (line[0] == '>' and str(id) == get_id(line)):
            return read_FASTA_sequence(file)


def search_FASTA_file_by_id(id, filename):
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
    basecomplement = {'A': 'T', 'C':'G', 'G':'C', 'T':'A'}
    letters = list(s)
    # find the complement
    letters = [basecomplement[base] for base in letters]
    
    return ''.join(letters)[::-1]

if __name__ == '__main__':
    # import pdb; pdb.set_trace() 
    id = 'chr01'
    seq = search_FASTA_file_by_id(id, 'yeast_chrom.fasta')
    #seq = search_FASTA_file_by_id(id, 'foo.txt')
    seq_frag = eval_Seq(seq, start = 230000, end = 230010, strand = '-')
    print(seq_frag)

    

