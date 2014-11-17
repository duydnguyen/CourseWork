1. Test 1
=========

python doViterbi.py Tests/transition1.txt Tests/emission1.txt 0 5 GCTT

output
Viterbi for state 1 time 1 : 0.1 maxstate 0
Viterbi for state 2 time 1 : 0.1 maxstate 0
Viterbi for state 3 time 1 : 0.0 maxstate 1
Viterbi for state 4 time 1 : 0.0 maxstate 2
Viterbi for state 1 time 2 : 0.004 maxstate 1
Viterbi for state 2 time 2 : 0.016 maxstate 2
Viterbi for state 3 time 2 : 0.018 maxstate 1
Viterbi for state 4 time 2 : 0.008 maxstate 2
Viterbi for state 1 time 3 : 0.0008 maxstate 1
Viterbi for state 2 time 3 : 0.00384 maxstate 2
Viterbi for state 3 time 3 : 0.0018 maxstate 3
Viterbi for state 4 time 3 : 0.00032 maxstate 2
Viterbi for state 1 time 4 : 0.00016 maxstate 1
Viterbi for state 2 time 4 : 0.0009216 maxstate 2
Viterbi for state 3 time 4 : 0.00018 maxstate 3
Viterbi for state 4 time 4 : 7.68e-05 maxstate 2
Viterbi probability: 9e-05
Viterbi path:
0
1 -> G
3 -> C
3 -> T
3 -> T


2. Test 2 
=========

python doViterbi.py Tests/transition2.txt Tests/emission2.txt 0 5 TAG

output
Viterbi for state 1 time 1 : 0.15 maxstate 0
Viterbi for state 2 time 1 : 0.2 maxstate 0
Viterbi for state 3 time 1 : 0.0 maxstate 1
Viterbi for state 4 time 1 : 0.0 maxstate 2
Viterbi for state 1 time 2 : 0.012 maxstate 1
Viterbi for state 2 time 2 : 0.064 maxstate 2
Viterbi for state 3 time 2 : 0.024 maxstate 1
Viterbi for state 4 time 2 : 0.004 maxstate 2
Viterbi for state 1 time 3 : 0.00048 maxstate 1
Viterbi for state 2 time 3 : 0.00512 maxstate 2
Viterbi for state 3 time 3 : 0.00288 maxstate 1
Viterbi for state 4 time 3 : 0.00512 maxstate 2
Viterbi probability: 0.004608
Viterbi path:
0
2 -> T
2 -> A
4 -> G

3. Test 3
=========

python doViterbi.py Tests/transition3.txt Tests/emission3.txt 0 3 GCTT

output
Viterbi for state 1 time 1 : 0.1 maxstate 0
Viterbi for state 2 time 1 : 0.0 maxstate 1
Viterbi for state 1 time 2 : 0.008 maxstate 1
Viterbi for state 2 time 2 : 0.008 maxstate 1
Viterbi for state 1 time 3 : 0.00256 maxstate 1
Viterbi for state 2 time 3 : 0.00016 maxstate 1
Viterbi for state 1 time 4 : 0.0008192 maxstate 1
Viterbi for state 2 time 4 : 5.12e-05 maxstate 1
Viterbi probability: 4.608e-05
Viterbi path:
0
1 -> G
1 -> C
1 -> T
2 -> T
