LearnMotif: Implemetation of MEME OOPS algorithm.
==================================================

1.  LearnMotif takes as input a set of DNA sequences and a width W, 
and learns an OOPS motif model for a motif of width W.

2. To run the program

+ The command line to run LearnMotif should be

'python LearnMotif.py sequences file width model file positions file num_seed max_iteration epsilon'


++ ***sequences_file*** is the file containing the input sequences. This file contains the DNA sequences, with one sequence per line.

++ ***width*** is the width of the motif model to learn

++ ***model_file*** is the name of a file to which you will output the learned motif model. This file contains a tab-delimited profile matrix (PWD matrix), with the background frequencies in the first
column.

++ ***positions_file*** is the name of a file to which you will output the predicted location of the motif in each sequences. This file simply contains a list of the best position for
the motif in each sequence, one position per line. 

++ ***num_seed*** is the maximum number of seeds allowed. Each value of seed corresponds to a different starting point (i.e., different PWD matrix)

++ ***max_iteration*** is the maximum number of iterations allowed in the E-M algorithm.

++ ***epsilon*** is the cut-off for change in Log Likelihood
 
+ To run question 2, the sequence_file = ***hw1_hidden_motif.txt*** should be in the same path as file ***LearnMotif.py***. For ***num_seed*** = 100, ***max_iteration*** = 10000, ***epsilon*** = 0.001, simply run the following command line

`python LearnMotif.py hw1_hidden_motif.txt 14 model_file positions_file 10 10000 0.001'