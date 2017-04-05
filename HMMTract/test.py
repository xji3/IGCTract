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
    for paralog in pairs[:]:
        print
        print paralog
        alignment_file = './MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
        MG94_IGC = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
        IGC_sitewise_lnL_file = './summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
        Force_sitewise_lnL_file = './summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
        Total_blen = sum([MG94_IGC.edge_to_blen[edge] for edge in MG94_IGC.edge_list if edge != ('N0', 'kluyveri')])

        test = HMMTract(IGC_sitewise_lnL_file, Force_sitewise_lnL_file, state_list,
                        Total_blen, MG94_IGC.tau)

        self = test

        #x = np.log([1.9, 0.3]) 
        #lnL_array = test.Forward(True, x)
        result = test.get_mle(False)

    #    test.IGC_sitewise_lnL = np.array(test.IGC_sitewise_lnL) * 0.001

        lnL_array, state_array = test.Viterbi()
        print 'Estimated average tract length: ', 3.0 / self.tract_p
        print 'Distn of S: ', test.get_marginal_state_distn()
        #print 'lnPr(x) - lnPr(x|S=0)', test.IGC_sitewise_lnL - test.Emi[0, :]
        print 'Viterbi path: ', state_array[list(lnL_array[:, -1]).index(max(lnL_array[:, -1]))]
    
