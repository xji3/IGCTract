from HMMTract import *
from IGCexpansion.CodonGeneconv import ReCodonGeneconv

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

        test = HMMTract(IGC_sitewise_lnL_file, Force_sitewise_lnL_file, state_list,
                        Total_blen, MG94_IGC.tau)
        MG94_IGC_lnL = -test.objective_1D(False, [0.0])

        print 
        print 'Tract length 1 lnL: ', MG94_IGC_lnL, output_ctrl
        result = test.get_mle(False)
        summary_mat.append([MG94_IGC_lnL, -result['fun'], 3.0 / test.tract_p,
                            test.get_marginal_state_distn()[0], test.get_marginal_state_distn()[1]])

        lnL_array, state_array = test.Viterbi()
        print 'Maximum lnL: ', -result['fun'], output_ctrl
        print 'Estimated average tract length: ', 3.0 / test.tract_p, output_ctrl
        print 'Distn of S: ', test.get_marginal_state_distn(), output_ctrl
        print 'Viterbi path: ', ''.join([ str(item) for item in state_array[list(lnL_array[:, -1]).index(max(lnL_array[:, -1]))]]), output_ctrl

        lnL_arr = test.get_posterior()
        np.savetxt('./summary/' + '_'.join(paralog) + '_MG94_nonclock_HMM_log_posterior_ratio.txt', lnL_arr[1, :] - lnL_arr[0, :])
        np.savetxt('./summary/' + '_'.join(paralog) + '_MG94_nonclock_HMM_Viterbi_path.txt', state_array[list(lnL_array[:, -1]).index(max(lnL_array[:, -1]))])

        

    np.savetxt('./HMM_tract_MG94_nonclock_summary.txt', np.matrix(summary_mat), delimiter = '\t')
