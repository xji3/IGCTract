# This script is used to infer IGC tract using HMM
# Xiang Ji
# xji3@ncsu.edu

import numpy as np
import os

class HMMTract:
    def __init__(self, IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                 State_List):
        self.IGC_sitewise_lnL   = self.read_lnL(IGC_sitewise_lnL_file)
        self.Force_sitewise_lnL = self.read_lnL(Force_sitewise_lnL_file)
        self.StateList          = State_List

        # Now inferrence related parameters
        self.Ptr = None        # Transition probability matrix
        self.Emi = None        # Emission probability matrix
        self.eta = None        # IGC initiation rate
        self.tract_p = None    # IGC tract distribution p as in Geo(p)


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
                
        
