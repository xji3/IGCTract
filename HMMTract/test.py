from HMMTract import *
from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import numdifftools as nd

if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))


    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    newicktree = './YeastTree.newick'
    Force = None
    output_ctrl = ''
    summary_mat = []
    for paralog in pairs[:]:
        print 
        print '**' + '_'.join(paralog)+ '**', output_ctrl
        
        alignment_file = './MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
        MG94_IGC = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
        IGC_sitewise_lnL_file = './summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
        Force_sitewise_lnL_file = './summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
        Total_blen = sum([MG94_IGC.edge_to_blen[edge] for edge in MG94_IGC.edge_list if edge != ('N0', 'kluyveri')])
        seq_index_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_seq_index.txt'

        test = HMMTract(IGC_sitewise_lnL_file, Force_sitewise_lnL_file, state_list,
                        Total_blen, MG94_IGC.tau, seq_index_file)
        
        
        MG94_IGC_lnL = -test.objective_1D(False, [0.0])

        print 
        print 'Tract length 1 lnL: ', MG94_IGC_lnL, output_ctrl
        if not "YNL069C" in paralog:
            test.objective_1D(False, [np.log(3.0 / 200.0)])
        result = test.get_mle(False)

        # Now use numdifftools to get Hessian (it's rather 2nd derivative though)
        #f = nd.Hessian(partial(test.objective_1D, False))
        f  = nd.Derivative(partial(test.objective_1D, False), n = 1)
        ff = nd.Derivative(partial(test.objective_1D, False), n = 2)
        #fisher_info = f(test.x[1:])[0,0]
        
        first_deriv  = -f(test.x[1:])
        second_deriv = -ff(test.x[1:])

        # store summary values
        summary_mat.append([MG94_IGC_lnL, -result['fun'], 3.0 / test.tract_p,
                            test.get_marginal_state_distn()[0], test.get_marginal_state_distn()[1],
                            first_deriv, second_deriv])
        
        lnL_array, Viterbi_path = test.Viterbi()
        print 'Maximum lnL: ', -result['fun'], output_ctrl
        print 'Estimated average tract length: ', 3.0 / test.tract_p, output_ctrl
        print 'Distn of S: ', test.get_marginal_state_distn(), output_ctrl
        print 'Viterbi path: ', ''.join([str(item) for item in Viterbi_path]), output_ctrl

        lnL_arr = test.get_posterior()
        np.savetxt('./summary/' + '_'.join(paralog) + '_MG94_nonclock_HMM_log_posterior_ratio.txt', lnL_arr[1, :] - lnL_arr[0, :])
        np.savetxt('./summary/' + '_'.join(paralog) + '_MG94_nonclock_HMM_Viterbi_path.txt', Viterbi_path)

        if True:
            lnL_surf = []
            tract_p_list = np.log(3.0 / np.arange(3, 501))
            for ln_tract_p in tract_p_list:
                lnL_surf.append(test.objective_1D(False, [ln_tract_p]))
            np.savetxt('./summary/' + '_'.join(paralog) + '_MG94_nonclock_HMM_lnL_surface.txt', np.array(lnL_surf))


    np.savetxt('./HMM_tract_MG94_nonclock_summary.txt', np.matrix(summary_mat), delimiter = '\t')
