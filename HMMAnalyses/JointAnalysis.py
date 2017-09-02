# This file is for joint analysis of MG94+IS-IGC+HMM over all gene pairs
# such that they share one common tract length distribution parameter tract_p
# It's kept out of the IGCexpansion package.
# xji3@ncsu.edu
from IGCexpansion.HMMTract import *
from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
import numdifftools as nd
from copy import deepcopy

class JointHMM:
    auto_save_step = 3
    def __init__(self, save_file,
                 # IndCodonGeneconv inputs
                 newicktree, alignment_file_list, paralog_list,
                 summary_path, save_path,
                 # parameter storage
                 x,
                 # HMMTract inputs
                 IGC_sitewise_lnL_file_list, NOIGC_sitewise_lnL_file_list,
                 State_list, seq_index_file_list,
                 model = 'MG94',
                 force = None, nsites = None, clock = None):
        # x is an array of length 2
        # x = np.log([tract_p])
        self.x = x
        
        # Store some info
        self.IGC_sitewise_lnL_file_list = IGC_sitewise_lnL_file_list
        self.NOIGC_sitewise_lnL_file_list = NOIGC_sitewise_lnL_file_list
        self.seq_index_file_list = seq_index_file_list

        # Now check length
        assert(len(alignment_file_list) == len(paralog_list) == len(IGC_sitewise_lnL_file_list) \
               == len(NOIGC_sitewise_lnL_file_list) == len(seq_index_file_list))
        self.MG94_IGC_list = [IndCodonGeneconv(newicktree, alignment_file_list[i], paralog_list[i], Model = model, \
                                               Force = force, clock = clock, \
                                               save_name = save_path + '_'.join(paralog_list[i]) + '_MG94_IGC_JointHMM_save.txt')\
                              for i in range(len(alignment_file_list))]
        self.hmmtract_list = list()
        self.save_path = save_path
        self.save_file = save_file

        self.auto_save_step = 0
        self.initialize()

    def initialize(self):
        for i in range(len(self.MG94_IGC_list)):
            MG94_IGC = self.MG94_IGC_list[i]
            MG94_IGC.get_mle(True, True, 0, 'BFGS')
            MG94_IGC.save_x()
            MG94_IGC.get_sitewise_loglikelihood_summary(self.IGC_sitewise_lnL_file_list[i], False)
            MG94_IGC.get_sitewise_loglikelihood_summary(self.NOIGC_sitewise_lnL_file_list[i], True)
            outgroup_branch = [edge for edge in MG94_IGC.edge_list if edge[0] == 'N0' and edge[1] != 'N1'][0]
            Total_blen = sum([MG94_IGC.edge_to_blen[edge] for edge in MG94_IGC.edge_list if edge != outgroup_branch])
            hmmtract = HMMTract(IGC_sitewise_lnL_file_list[i], NOIGC_sitewise_lnL_file_list[i], State_list, Total_blen,
                                MG94_IGC.tau, seq_index_file_list[i])
            self.hmmtract_list.append(hmmtract)

    def _loglikelihood(self, display, x):
        assert(len(x) == 1)
        self.x = x
        ll = 0.0
        for i in range(len(self.MG94_IGC_list)):
            hmmtract = self.hmmtract_list[i]
            ll += -hmmtract.objective_1D(False, self.x)

        if display:
            print(ll, np.exp(self.x[0]), self.x[0])
        self.auto_save_step += 1
        return ll

    def objective(self, display, x):
        return -self._loglikelihood(display, x)
                              
    def save_x(self):
        np.savetxt(open(self.save_file, 'w+'), self.x.T)
        for MG94_IGC in self.MG94_IGC_list:
            MG94_IGC.save_x()

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
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))


    
    newicktree = './YeastTree.newick'
    State_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    paralog_list = pairs
    save_file = './save/' + str(len(paralog_list)) + '_pair_joint_HMM_save.txt'
    alignment_file_list = ['../HMMTract/MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'\
                           for paralog in paralog_list]
    
    summary_path = './summary/'
    save_path = './save/'
    x = [0.0]
    IGC_sitewise_lnL_file_list = ['./summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt' \
                                  for paralog in paralog_list]
    NOIGC_sitewise_lnL_file_list = ['./summary/NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt' \
                                    for paralog in paralog_list]
    seq_index_file_list = ['../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_seq_index.txt' \
                           for paralog in paralog_list]

    test = JointHMM(save_file, newicktree, alignment_file_list, paralog_list, summary_path, save_path,
                    x, IGC_sitewise_lnL_file_list, NOIGC_sitewise_lnL_file_list, State_list, seq_index_file_list)
    self = test
    print (test._loglikelihood(True, test.x))
    test.get_mle()
    print (test.get_Hessian())

    log_p_list = np.log(3.0/np.array(range(3, 1001)))
    plot_file = './plot/' + str(len(paralog_list)) + '_pair_joint_HMM_lnL_1D_surface.txt'
    test.plot_tract_p(log_p_list, plot_file)

    
