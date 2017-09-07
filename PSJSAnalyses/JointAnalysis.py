# This file is for joint analysis of MG94+IS-IGC+HMM over all gene pairs
# such that they share one common tract length distribution parameter tract_p
# It's kept out of the IGCexpansion package.
# xji3@ncsu.edu
from IGCexpansion.PSJSGeneconv import *
from IGCexpansion.JSGeneconv import JSGeneconv
import argparse
from collections import namedtuple
import numpy as np
import numdifftools as nd
from copy import deepcopy

class JointPSJS:
    auto_save_step = 2
    def __init__(self, save_file,
                 # JSGeneconv inputs
                 newicktree, alignment_file_list, gene_to_orlg_file_list,
                 DupLosList, paralog_list, 
                 rate_variation, node_to_pos, terminal_node_list, JS_save_file_list,
                 # parameter storage
                 x,
                 # PSJSGeneconv additional inputs
                 IGC_pm, seq_index_file_list, PSJS_save_file_list, PSJS_log_file_list,
                 model = 'HKY',
                 force = None, nsites = None, clock = None):
        # x is an array of length 1
        # x = np.log([tract_p])
        self.x = x
        
        # Store some info
        self.seq_index_file_list = seq_index_file_list

        
        if rate_variation:
            x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
        else:
            x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])

        # Now check length
        assert(len(alignment_file_list) == len(paralog_list) == len(gene_to_orlg_file_list) \
               == len(JS_save_file_list) == len(seq_index_file_list) == len(PSJS_save_file_list) == len(PSJS_log_file_list))

        # Now assign list of class instances
        self.JS_test_list = [JSGeneconv(alignment_file_list[i], gene_to_orlg_file_list[i], True, newicktree,
                                        DupLosList, x_js, model, IGC_pm, rate_variation, node_to_pos,
                                        terminal_node_list, JS_save_file_list[i])\
                             for i in range(len(alignment_file_list))]

        x_psjs = np.concatenate((x_js[:-1], [x_js[-1] - np.log(100.0), - np.log(100.0) ]))
        self.PSJS_test_list = [PSJSGeneconv(alignment_file_list[i], gene_to_orlg_file_list[i], seq_index_file_list[i],
                                            True, True, newicktree, DupLosList, x_psjs, model, IGC_pm,
                                            rate_variation, node_to_pos, terminal_node_list, PSJS_save_file_list[i],
                                            PSJS_log_file_list[i])\
                               for i in range(len(alignment_file_list))]

        self.save_file = save_file

        self.auto_save_step = 0
        self.initialize()

    def initialize(self):
        for i in range(len(self.JS_test_list)):
            test_JS = self.JS_test_list[i]
            test_JS.get_mle()
            test_JS.save_x()


            PSJS_IGC = self.PSJS_test_list[i]
            x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(100.0), - np.log(100.0) ],
                           test_JS.x[len(test_JS.jsmodel.x_js):]))
            PSJS_IGC.unpack_x(x)
            

    def _loglikelihood(self, display, x):
        assert(len(x) == 1)
        self.x = x
        ll = 0.0
        for i in range(len(self.PSJS_test_list)):
            PSJS_IGC = self.PSJS_test_list[i]
            ll += -PSJS_IGC.objective_tract_p(False, self.x[0])

        if display:
            print(ll, np.exp(self.x[0]), self.x[0])
        self.auto_save_step += 1
        return ll

    def objective(self, display, x):
        return -self._loglikelihood(display, x)
                              
    def save_x(self):
        np.savetxt(open(self.save_file, 'w+'), self.x.T)
        for i in range(len(self.JS_test_list)):
            test_JS = self.JS_test_list[i]
            test_JS.save_x()

            PSJS_IGC = self.PSJS_test_list[i]
            PSJS_IGC.save_x()

    def get_mle(self, display = True):
        f = partial(self.objective, display)
        bnds = [(None, 0.0)]
        guess_x = deepcopy(self.x)
        result = scipy.optimize.minimize(f, guess_x, jac = False, method = 'L-BFGS-B', bounds = bnds)
        self.save_x()

        if display:
            print(result)

        return result


    def plot_tract_p(self, log_p_list, plot_file):
        ll_list = []
        for log_p in log_p_list:
            ll = self._loglikelihood(False, [log_p])
            ll_list.append(ll)

        with open(plot_file, 'w+') as f:
            f.write('# log_p \t lnL \t \n')
            for it in range(len(log_p_list)):
                f.write('\t'.join([str(log_p_list[it]), str(ll_list[it]), '\n']))

    def get_Hessian(self):
        f = nd.Derivative(partial(self._loglikelihood, False), n = 1)
        ff = nd.Derivative(partial(self._loglikelihood, False), n = 2)
        grad = -f(self.x)
        hess = -ff(self.x)
        return grad, hess
    
if __name__ == '__main__':
    pairs = []
    all_pairs = '../Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))


    
    newicktree = '../YeastTree.newick'
    
    paralog_list = pairs
    
    save_file = './save/' + str(len(paralog_list)) + '_pair_joint_PSJS_save.txt'
    alignment_file_list = ['../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'\
                           for paralog in paralog_list]
    gene_to_orlg_file_list = ['../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt' for paralog in paralog_list]
    DupLosList = '../YeastTestDupLost.txt'
    rate_variation = True
    IGC_pm = 'One rate'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    if rate_variation:
        JS_save_file_list = ['./save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_JointAnalysis_save.txt' \
                             for paralog in paralog_list]
        PSJS_save_file_list = ['./save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_JointAnalysis_save.txt' \
                             for paralog in paralog_list]
        PSJS_log_file_list = ['./log/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_JointAnalysis_log.txt' \
                             for paralog in paralog_list]
    else:
        JS_save_file_list = ['./save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_JointAnalysis_save.txt' \
                             for paralog in paralog_list]
        PSJS_save_file_list = ['./save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_JointAnalysis_save.txt' \
                             for paralog in paralog_list]
        PSJS_log_file_list = ['./log/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_JointAnalysis_log.txt' \
                             for paralog in paralog_list]
        
    x = [0.0]
    
    seq_index_file_list = ['../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_seq_index.txt' \
                           for paralog in paralog_list]

    test = JointPSJS(save_file,
                 newicktree, alignment_file_list, gene_to_orlg_file_list,
                 DupLosList, paralog_list, 
                 rate_variation, node_to_pos, terminal_node_list, JS_save_file_list,
                 # parameter storage
                 x,
                 # PSJSGeneconv additional inputs
                 IGC_pm, seq_index_file_list, PSJS_save_file_list, PSJS_log_file_list,)
    self = test
#    print (test._loglikelihood(True, test.x))
    test.get_mle()
##    print (test.get_Hessian())
##
##    log_p_list = np.log(3.0/np.array(range(3, 1001)))
##    plot_file = './plot/' + str(len(paralog_list)) + '_pair_joint_HMM_lnL_1D_surface.txt'
##    test.plot_tract_p(log_p_list, plot_file)

    
