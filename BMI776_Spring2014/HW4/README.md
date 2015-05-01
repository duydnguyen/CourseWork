#  Implemetations of Nussinov Algorithm

## Nussinov

1. *Nussinov* takes as input an RNA sequence, and output the base-paired positions of a structure that maximizes the number of base pairings.

2. To run the program

* The command line to run *Nussinov* should be

``python Nussinov.py RNASEQUENCE``

* RNASEQUENCE is a string of A, C, G, U characters. 

* The program prints to standard output the positions of all pairs of base-paired positions (using one-based indexing), in lexicographical order.

* For example, the command 

``python Nussinov.py UCCAGG``

should output

2 6 

3 5
