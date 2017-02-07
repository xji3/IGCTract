from IGCexpansion.PSJSGeneconv import *
from IGCexpansion.JSGeneconv import *
import argparse
from collections import namedtuple
import numpy as np

def main(args):
    paralog = [args.paralog1, args.paralog2]
    
    gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt'
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_input.fasta'

    tree_newick = '../YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])
    IGC_pm = 'One rate'

    save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
    summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'

    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
                      node_to_pos, terminal_node_list, save_file)
    test_JS.get_mle()
    test_JS.get_expectedNumGeneconv()
    test_JS.get_individual_summary(summary_file)

    if args.dim == 1:     
        save_file = './save/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
        summary_file = './summary/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
    elif args.dim == 2:
        save_file = './save/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
        summary_file = './summary/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
    else:
        save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
        summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
        
    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'

    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ] ))
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
                      node_to_pos, terminal_node_list, save_file)

    x = np.concatenate((x_js, test_JS.x[len(test_JS.jsmodel.x_js):]))
    test.unpack_x(x)
    if args.dim == 1:
        test.optimize_x_IGC(dimension = 1)
    elif args.dim == 2:
        test.optimize_x_IGC(dimension = 2)
    else:
        test.get_mle()
    test.get_individual_summary(summary_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--L', type = float, dest = 'tract_length', default = 30.0, help = 'Initial guess tract length')
    parser.add_argument('--D', type = int, dest = 'dim', default = 0, help = 'Dimension used in search with default value 0')
    
    main(parser.parse_args())

  

##    MyStruct = namedtuple('MyStruct', 'paralog1 paralog2 tract_length dim')
##    args = MyStruct(paralog1 = 'YLR406C', paralog2 = 'YDL075W', tract_length = 30.0, dim = 1)
##
##
##    paralog = [args.paralog1, args.paralog2]
##    
##    gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt'
##    alignment_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_input.fasta'
##
##    tree_newick = '../YeastTree.newick'
##    DupLosList = '../YeastTestDupLost.txt'
##    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
##    node_to_pos = {'D1':0}
##    pm_model = 'HKY'
##    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])
##    IGC_pm = 'One rate'
##
##    save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
##    summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
##
##    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
##                      node_to_pos, terminal_node_list, save_file)
##    test_JS.get_mle()
##    test_JS.get_expectedNumGeneconv()
##    test_JS.get_individual_summary(summary_file)
##
##    if args.dim == 1:     
##        save_file = './save/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
##        summary_file = './summary/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
##    elif args.dim == 2:
##        save_file = './save/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
##        summary_file = './summary/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
##    else:
##        save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
##        summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
##        
##    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'
##
##    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
##                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ] ))
##    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
##                      node_to_pos, terminal_node_list, save_file)
##
##    x = np.concatenate((x_js, test_JS.x[len(test_JS.jsmodel.x_js):]))
##    test.unpack_x(x)
##    if args.dim == 1:
##        test.optimize_x_IGC(dimension = 1)
##    elif args.dim == 2:
##        test.optimize_x_IGC(dimension = 2)
##    else:
##        test.get_mle()
##    test.get_individual_summary(summary_file)


##    pairs = []
##    all_pairs = './Finished_Pairs.txt'
##    with open(all_pairs, 'r') as f:
##        for line in f.readlines():
##            pairs.append(line.replace('\n','').split('_'))
##
##    for pair in pairs:  
