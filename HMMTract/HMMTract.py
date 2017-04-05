# This script is used to infer IGC tract using HMM
# Xiang Ji
# xji3@ncsu.edu

import numpy as np
import scipy, scipy.optimize, scipy.linalg 
from scipy.linalg import expm
from functools import partial
import os, sys

class HMMTract:
    def __init__(self, IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                 State_List, Total_blen, tau):
        self.IGC_sitewise_lnL   = self.read_lnL(IGC_sitewise_lnL_file)
        self.Force_sitewise_lnL = self.read_lnL(Force_sitewise_lnL_file)
        self.StateList          = State_List
        self.L                  = Total_blen

        # Now inferrence related parameters
        self.tau = tau         # estimated Tau value from MG94+IGC independent site model
        self.Ptr = None        # Transition probability matrix
        self.Emi = None        # Emission probability matrix
        self.eta = None        # IGC initiation rate
        self.tract_p = None    # IGC tract distribution p as in Geo(p)

        self.x      = None        # log array to store eta and tract_p values
        self.is_mle = False

        self.init_parameters()


    def init_parameters(self):
        self.update_by_x(np.log([self.tau * 0.05, 0.05]))
        assert(len(self.IGC_sitewise_lnL) == len(self.Force_sitewise_lnL))

    def read_lnL(self, sitewise_lnL_file):
        assert(os.path.isfile(sitewise_lnL_file))
        pos = []
        ll  = []
        with open(sitewise_lnL_file, 'rb') as f:
            for line in f:
                if line[0] == '#':
                    continue
                items = line.replace('\n', '').split('\t')
                pos.append(int(items[0]))
                ll.append(float(items[1]))

        return ll

    def update_by_x(self, x):
        assert(len(x) == 2)
        self.x = x
        self.eta, self.tract_p = np.exp(x)
        self.get_Ptr_analytical()
        self.get_Emi()
        
    def get_marginal_state_distn(self):
        assert( self.eta != None and self.tract_p != None)
        P_S0 = np.exp( -2.0 * self.eta / self.tract_p * self.L)
        P_S1 = 1.0 - P_S0
        return [P_S0, P_S1]

    def get_Ptr(self):
        p = self.tract_p
        Q = self.eta * np.matrix([[-1.0 - 1.0/p, 1.0,   1.0, 1.0/p - 1.0],
                                  [0.0, -1.0/p, 0.0, 1.0/p],
                                  [0.0, 0.0, -1.0/p, 1.0/p],
                                  [0.0, 0.0, 0.0, 0.0]], dtype = float)
        P_mat = scipy.linalg.expm2(2 * Q * self.L)
        
        distn = self.get_marginal_state_distn()
        Ptr = np.matrix([[P_mat[0, 0] / distn[0], P_mat[0, 1] / distn[0] ],
                         [P_mat[0, 2] / distn[1], P_mat[0, 3] / distn[1] ]])

        self.Ptr = np.log(Ptr)

    def get_Ptr_analytical(self):
        etl = np.exp(-2 * self.L * self.eta / self.tract_p)
        eel = np.exp(-2 * self.L * self.eta)
        Ptr = np.matrix([[eel, 1.0 - eel], [etl*(1.0 - eel)/(1.0 - etl), (1.0 - etl - (1.0 - eel) * etl)/(1.0 - etl)]], dtype = float)
        self.Ptr = np.log(Ptr)

    def get_Emi(self):
        self.Emi = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL)), dtype = float)
        distn = self.get_marginal_state_distn()
        for i in range(len(self.IGC_sitewise_lnL)):
            emission_0 = self.Force_sitewise_lnL[i]
            emission_1 = self.IGC_sitewise_lnL[i] + np.log(1.0 - np.exp(emission_0 - self.IGC_sitewise_lnL[i]) * distn[0]) - np.log(distn[1])

            if 1.0 - np.exp(emission_0 - self.IGC_sitewise_lnL[i]) * distn[0] < 0:
                print self.x
                sys.exit('something is wrong')

            self.Emi[:, i] = np.array([emission_0, emission_1])
            
        

    def Forward(self, display, x): # lnL by forward algorithm
        # update parameters first
        self.update_by_x(x)

        distn = self.get_marginal_state_distn()

        # Now create a 2 by nsites array for the dynamic programing
        lnL_array = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL) + 1), dtype = float)

        # Now add in initial distribution
        lnL_array[:, 0] = np.log(distn)

        # Now do the forward step
        for i in range(len(self.IGC_sitewise_lnL)):
            emission_0 = self.Emi[0, i]
            emission_1 = self.Emi[1, i]
            
            new_cond_lnL_0 = emission_0 + (lnL_array[0, i] + self.Ptr[0, 0]) + np.log(sum(np.exp([0.0, lnL_array[1, i] + self.Ptr[1, 0] - (lnL_array[0, i] + self.Ptr[0, 0])])))
            new_cond_lnL_1 = emission_1 + (lnL_array[0, i] + self.Ptr[0, 1]) + np.log(sum(np.exp([0.0, lnL_array[1, i] + self.Ptr[1, 1] - (lnL_array[0, i] + self.Ptr[0, 1])])))

            lnL_array[:, i + 1] = np.array([new_cond_lnL_0, new_cond_lnL_1])

        ll = sum(lnL_array[:, -1])
        if display:
            print '\t'.join([ str(item) for item in [ll] + list(np.exp(self.x))])

        return -ll

    def objective_1D(self, display, ln_p):
        x = np.array([np.log(self.tau) + ln_p[0], ln_p[0]])
        return self.Forward(display, x)
    

    def get_mle(self, display = True, OneDimension = True):
        self.update_by_x(self.x)
        if not OneDimension:
            f = partial(self.Forward, display)
            guess_x = self.x
            bnds = [(None, None), (None, 0.0)]

        else:
            f = partial(self.objective_1D, display)
            guess_x = self.x[1]
            bnds = [(None, 0.0)]

        result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)

        if display:
            print(result)
        if result['success']:
            self.update_by_x(self.x)
            self.is_mle = True
        return result

    def Viterbi(self):
        assert(self.is_mle)
        self.update_by_x(self.x)
        # same as the setup in Forward algorithm
        distn = self.get_marginal_state_distn()

        # Now create a 2 by nsites array for the Viterbi algorithm
        lnL_array = np.zeros((len(self.StateList), len(self.IGC_sitewise_lnL) + 1), dtype = float)
        state_array = [[], []]

        # Now add in initial distribution
        lnL_array[:, 0] = np.log(distn)

        # Now do the Viterbi algorithm
        for i in range(len(self.IGC_sitewise_lnL)):
            emission_0 = self.Emi[0, i]
            emission_1 = self.Emi[1, i]
            
            new_cond_lnL_0_list = [emission_0 + lnL_array[0, i] + self.Ptr[0, 0], emission_0 + lnL_array[1, i] + self.Ptr[1, 0]]
            lnL_array[0, i + 1] = max(new_cond_lnL_0_list)
            state_array[0].append(new_cond_lnL_0_list.index(lnL_array[0, i + 1]))

            new_cond_lnL_1_list = [emission_1 + lnL_array[0, i] + self.Ptr[0, 1], emission_1 + lnL_array[1, i] + self.Ptr[1, 1]]
            lnL_array[1, i + 1] = max(new_cond_lnL_1_list)
            state_array[1].append(new_cond_lnL_1_list.index(lnL_array[1, i + 1]))



##        if display:
##            print

        return lnL_array, state_array
        
        
            
        

                
        
