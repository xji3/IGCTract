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
    
    IGC_pm = 'One rate'

    if args.rate_variation:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    else:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])


    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, args.cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                         args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
    test_JS.get_mle()
    #test_JS.get_expectedNumGeneconv()
    test_JS.get_individual_summary(summary_file)

    if args.dim == 1:     
        save_file = './save/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
        log_file = './log/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_log.txt'
        summary_file = './summary/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
    elif args.dim == 2:
        save_file = './save/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
        log_file = './log/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_log.txt'
        summary_file = './summary/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
    else:
        save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
        log_file = './log/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_log.txt'
        summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'

    if args.rate_variation:
        if args.allow_same_codon:
            save_file = save_file.replace('_nonclock', '_rv_SCOK_nonclock')
            log_file = log_file.replace('_nonclock', '_rv_SCOK_nonclock')
            summary_file = summary_file.replace('_nonclock', '_rv_SCOK_nonclock')
        else:
            save_file = save_file.replace('_nonclock', '_rv_NOSC_nonclock')
            log_file = log_file.replace('_nonclock', '_rv_NOSC_nonclock')
            summary_file = summary_file.replace('_nonclock', '_rv_NOSC_nonclock')


        
    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'

    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ]))
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, args.cdna, args.allow_same_codon, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                      args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
##    x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
##                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ],
##                           test_JS.x[len(test_JS.jsmodel.x_js):]))
##    test.unpack_x(x)

    if args.dim == 1:
        test.optimize_x_IGC(dimension = 1)
    elif args.dim == 2:
        test.optimize_x_IGC(dimension = 2)
    else:
        test.get_mle()
    test.get_individual_summary(summary_file)

    pairwise_lnL_summary_file = summary_file.replace('_summary.txt', '_lnL_summary.txt')
    test.get_pairwise_loglikelihood_summary(pairwise_lnL_summary_file)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--L', type = float, dest = 'tract_length', default = 30.0, help = 'Initial guess tract length')
    parser.add_argument('--D', type = int, dest = 'dim', default = 0, help = 'Dimension used in search with default value 0')
    parser.add_argument('--heterogeneity', dest = 'rate_variation', action = 'store_true', help = 'rate heterogeneity control')
    parser.add_argument('--homogeneity', dest = 'rate_variation', action = 'store_false', help = 'rate heterogeneity control')
    parser.add_argument('--coding', dest = 'cdna', action = 'store_true', help = 'coding sequence control')
    parser.add_argument('--noncoding', dest = 'cdna', action = 'store_false', help = 'coding sequence control')
    parser.add_argument('--samecodon', dest = 'allow_same_codon', action = 'store_true', help = 'whether allow pair sites from same codon')
    parser.add_argument('--no-samecodon', dest = 'allow_same_codon', action = 'store_false', help = 'whether allow pair sites from same codon')
    
    main(parser.parse_args())

  

##    MyStruct = namedtuple('MyStruct', 'paralog1 paralog2 tract_length dim cdna rate_variation allow_same_codon')
##    args = MyStruct(paralog1 = 'YDR418W', paralog2 = 'YEL054C', tract_length = 30.0, dim = 1, cdna = True, rate_variation = False, allow_same_codon = False)
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
##    
##    IGC_pm = 'One rate'
##
##    if args.rate_variation:
##        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
##        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
##        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
##    else:
##        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
##        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
##        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])
##
##
##    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, args.cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
##                         args.rate_variation, node_to_pos, terminal_node_list, save_file)
##    test_JS.get_mle()
##    #test_JS.get_expectedNumGeneconv()
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
##    if args.rate_variation:
##        if args.allow_same_codon:
##            save_file = save_file.replace('_nonclock', '_rv_SCOK_nonclock')
##            summary_file = summary_file.replace('_nonclock', '_rv_SCOK_nonclock')
##        else:
##            save_file = save_file.replace('_nonclock', '_rv_NOSC_nonclock')
##            summary_file = summary_file.replace('_nonclock', '_rv_NOSC_nonclock')
##
##
##        
##    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'
##
##    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
##                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ]))
##    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, args.cdna, args.allow_same_codon, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
##                      args.rate_variation, node_to_pos, terminal_node_list, save_file)
##    x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
##                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ],
##                           test_JS.x[len(test_JS.jsmodel.x_js):]))
##    test.unpack_x(x)
##
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
