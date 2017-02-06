from IGCexpansion.PSJSGeneconv import *
import os, argparse
from collections import namedtuple
from copy import deepcopy
import numpy as np

if __name__ == '__main__':
    pairs = []
    all_pairs = './Finished_Pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    tract_length_list = np.log(np.arange(1.0, 20.0, 1.0))
    
    for pair in [('YLR406C', 'YDL075W')]:#('YLR406C', 'YDL075W')]:    

        MyStruct = namedtuple('MyStruct', 'paralog1 paralog2 tract_length')
        args = MyStruct(paralog1 = pair[0], paralog2 = pair[1], tract_length = 30.0)


        paralog = [args.paralog1, args.paralog2]
        
        gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt'
        alignment_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_input.fasta'

        tree_newick = '../YeastTree.newick'
        DupLosList = '../YeastTestDupLost.txt'
        terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
        node_to_pos = {'D1':0}
        pm_model = 'HKY'
        IGC_pm = 'One rate'

        save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'

        x_js = np.loadtxt(save_file)[:6]
        lnTau = x_js[-2] - x_js[-1]
        test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
                          node_to_pos, terminal_node_list, save_file)
        test.unpack_x(test.x)
        #test.get_mle()
        #test.get_individual_summary(summary_file)
        x_plot = deepcopy(test.x)
        if not os.path.isdir('./plot'):
            os.mkdir('./plot')
        plot_file_name = './plot/' + '_'.join(pair) + '_PSJS_lnL_TractLength.txt'
        #plot_file_name = './plot/' + '_'.join(pair) + '_PSJS_lnL_TractLength_zoomed.txt'
        with open(plot_file_name, 'w+') as f:
            tract_length = np.exp(-x_js[-1])
            #lnL = test.objective_wo_gradient(True, test.x)
            f.write('\t'.join(['#Tract_length', 'Init_rate', 'lnL', '\n']))

        for tract_length in tract_length_list:
            init_rate_list = lnTau - tract_length + np.log(np.arange(0.5, 2.0, 0.1))
            for init_rate in init_rate_list:
                x_plot[5] = -tract_length
                x_plot[4] = init_rate
                test.unpack_x(x_plot)
                lnL = test.objective_wo_gradient(True, x_plot)
                with open(plot_file_name, 'a') as f:
                    f.write('\t'.join([str(np.exp(tract_length)), str(np.exp(init_rate)), str(lnL), '\n']))
            
        
