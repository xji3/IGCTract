##This code is for ST590C HW2
##Xiang Ji
##xji3@ncsu.edu
from app import *

state_list = ['AT rich','GC rich']
inits = array([1.0, 0.0])
read_dir = './'


##################################################################################
################ Question 1 ######################################################
##################################################################################
P1 = mat([[0.98, 0.02],[0.05, 0.95]])
E = array([[0.3, 0.2, 0.2, 0.3],[0.2, 0.3, 0.3, 0.2]])
Q1 = HMM(P1, E, state_list, read_dir, read_dir)
logp1 = Q1.Forward(inits)
print 'The log-likelihood value of Question 1 is ', logp1
print

##################################################################################
################ Question 2 ######################################################
##################################################################################
P2 = mat([[0.8, 0.2],[0.5, 0.5]])
Q2 = HMM(P2, E, state_list, read_dir, read_dir)
logp2 = Q2.Forward(inits)
print 'The log-likelihood value of Question 2 is ', logp2
print 

##################################################################################
################ Question 3 ######################################################
##################################################################################
P3 = mat([[0.51, 0.49],[0.49, 0.51]])
Q3 = HMM(P3, E, state_list, read_dir, read_dir)
logp3 = Q3.Forward(inits)
print 'The log-likelihood value of Question 3 is ', logp3
print

##################################################################################
################ Question 4 ######################################################
##################################################################################
score4, path4 = Q1.Viterbi(inits)
print 'Best Viterbi log Score for question 4 is ',score4
Q1.printstatelist()
displaypath(path4)
print

##################################################################################
################ Question 5 ######################################################
##################################################################################
score5, path5 = Q2.Viterbi(inits)
print 'Best Viterbi log Score for question 5 is ',score5
Q2.printstatelist()
displaypath(path5)
print

##################################################################################
################ Question 6 ######################################################
##################################################################################
score6, path6 = Q3.Viterbi(inits)
print 'Best Viterbi log Score for question 6 is ',score6
Q3.printstatelist()
displaypath(path6)

test = True
for i in range(1,len(path6)):
    if Q3.seq[i]=='A' or Q3.seq[i]=='T':
            st = 0
    else:
            st = 1
    
    test = test and (path6[i]==st)
    if test == False:
            print i

if test:
    print 'The path is as expected'
else:
    print 'The path is not as expected'
