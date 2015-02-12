LearnMotif: Implemetation of MEME OOPS algorithm.
==================================================

1.  LearnMotif takes as input a set of DNA sequences and a width W, 
and learns an OOPS motif model for a motif of width W.

2. To run the program

+ The command line to run LearnMotif should be

`python LearnMotif.py sequences_file width model_file positions_file`

where ***sequences_file*** is the file containing the input sequences, ***width*** is the width
of the motif model to learn, *model_file* is the name of a file to which you will output the learned motif model, and ***positions_file*** is the name of a file to which you will output the predicted location of the motif in each sequences.

+ The input ***sequences_file*** will contain the DNA sequences, with one sequence per
line.

The output positions file should simply contain a list of the best position for
the motif in each sequence, one position per line. The ***model_file*** should contain a tab-delimited profile matrix (PWD matrix), with the background frequencies in the first
column.