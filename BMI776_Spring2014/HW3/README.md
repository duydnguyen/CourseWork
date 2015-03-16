 Implemetations of maximal dependence decomposition(MDD), interpolated Markov models (IMMs), and position weight matrix (PWD)
==================================================

For a better view of this file in Markdown, please go to `https://github.com/duydnguyen/CourseWork/tree/master/BMI776_Spring2014/HW3`

pwm_predict
-----------


1. *pwm_predict* learns PWM models from training data and output scores for test data.

2. To run the program

* The command line to run *pwm_predict* should be

``python pwm_predict.py train_real train_false test test.scores``

+ ***train_real*** is the file containing a set of positive train examples.

+ ***train_false*** is the file containing a set of negative train examples.

+ ***test*** is a file containing a set of test sequences for which you are to output scores.

+ ***test.scores*** is the file into which scores will be written.

* To run *pwm_predict* for the current Data set, simply run the following command line (assuming *pwm_predict.py* is in the same directory with data files)

``python pwm_predict.py hw3_train_real hw3_train_false hw3_test_real test_scores_real``

``python pwm_predict.py hw3_train_real hw3_train_false hw3_test_false test_scores_false``