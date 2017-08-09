from IGCexpansion.PSJSGeneconv import *
from IGCexpansion.JSGeneconv import *
import os, argparse
from collections import namedtuple
from copy import deepcopy
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
                         args.rate_variation, node_to_pos, terminal_node_list, save_file)
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

    lnTau = test_JS.jsmodel.x_js[-1]
    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ lnTau - np.log(args.tract_length), - np.log(args.tract_length) ]))
    
    
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, args.cdna, args.allow_same_codon, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                      args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
    x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ],
                           test_JS.x[len(test_JS.jsmodel.x_js):]))
    test.unpack_x(x)


    #test.get_mle()
    #test.get_individual_summary(summary_file)
    x_plot = deepcopy(test.x)
    if not os.path.isdir('./plot'):
        os.mkdir('./plot')

    if not os.path.isdir('./plot/' + '_'.join(paralog)):
        os.mkdir('./plot/' + '_'.join(paralog))
        
    if args.zoom_in:
        off_ratio = 0.1
    else:
        off_ratio = 0.5
        
    tract_length_start = np.exp(-x_js[-1]) * (1.0 - off_ratio)
    tract_length_end   = np.exp(-x_js[-1]) * (1.0 + off_ratio)
    tract_length_space = (tract_length_end - tract_length_start)/5
    tract_length_list = np.log(np.arange(tract_length_start, tract_length_end + tract_length_space, tract_length_space))

    plot_file_name = './plot/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_PSJS_HKY_lnL_'
    if args.rate_variation:
        if args.allow_same_codon:
            plot_file_name += 'rv_SCOK_nonclock_'
        else:
            plot_file_name += 'rv_NOSC_nonclock_'
    else:
        plot_file_name += 'nonclock_'

    plot_file_name += 'dim_' + str(args.dim) + '_'
            
    if args.zoom_in:
        plot_file_name += 'offratio_' + str(off_ratio) + '_zoomed_in.txt'
    else:
        plot_file_name += 'offratio_' + str(off_ratio) +'.txt'
    with open(plot_file_name, 'w+') as f:
        tract_length = np.exp(-x_js[-1])
        #lnL = test.objective_wo_gradient(True, test.x)
        f.write('\t'.join(['#Tract_length', 'Init_rate', 'lnL', '\n']))



    for tract_length in tract_length_list:
        init_rate_list = lnTau - tract_length + np.log(np.arange((1.0 - off_ratio), (1.0 + off_ratio), 2 * off_ratio / 4.0))
        for init_rate in init_rate_list:
            x_plot[5] = -tract_length
            x_plot[4] = init_rate
            test.unpack_x(x_plot)
            lnL = test.objective_wo_gradient(True, x_plot)
            with open(plot_file_name, 'a') as f:
                f.write('\t'.join([str(np.exp(tract_length)), str(np.exp(init_rate)), str(lnL), '\n']))


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
    parser.add_argument('--zoom', dest = 'zoom_in', action = 'store_true', help = 'clock control')
    parser.add_argument('--no-zoom', dest = 'zoom_in', action = 'store_false', help = 'clock control')
    
    main(parser.parse_args())        


##    MyStruct = namedtuple('MyStruct', 'paralog1 paralog2 tract_length dim cdna rate_variation allow_same_codon zoom_in')
##    args = MyStruct(paralog1 = 'YDR418W', paralog2 = 'YEL054C', tract_length = 30.0, dim = 1, cdna = True, rate_variation = True, allow_same_codon = True,
##                    zoom_in = True)
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
##        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
##        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
##        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
##    else:
##        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
##        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
##        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
##        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])
##
##    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, args.cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
##                         args.rate_variation, node_to_pos, terminal_node_list, save_file)
##    test_JS.get_mle()
##    #test_JS.get_expectedNumGeneconv()
##    test_JS.get_individual_summary(summary_file)
##
##    if args.dim == 1:     
##        save_file = './save/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
##        log_file = './log/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_log.txt'
##        summary_file = './summary/PSJS_dim_1_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
##    elif args.dim == 2:
##        save_file = './save/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
##        log_file = './log/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_log.txt'
##        summary_file = './summary/PSJS_dim_2_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
##    else:
##        save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_save.txt'
##        log_file = './log/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_log.txt'
##        summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_init_' + str(args.tract_length) + '_nonclock_summary.txt'
##
##    if args.rate_variation:
##        if args.allow_same_codon:
##            save_file = save_file.replace('_nonclock', '_rv_SCOK_nonclock')
##            log_file = log_file.replace('_nonclock', '_rv_SCOK_nonclock')
##            summary_file = summary_file.replace('_nonclock', '_rv_SCOK_nonclock')
##        else:
##            save_file = save_file.replace('_nonclock', '_rv_NOSC_nonclock')
##            log_file = log_file.replace('_nonclock', '_rv_NOSC_nonclock')
##            summary_file = summary_file.replace('_nonclock', '_rv_NOSC_nonclock')
##
##    
##        
##    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'
##
##    lnTau = test_JS.jsmodel.x_js[-1]
##    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
##                           [ lnTau - np.log(args.tract_length), - np.log(args.tract_length) ]))
##    
##    
##    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, args.cdna, args.allow_same_codon, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
##                      args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
##    x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
##                           [ test_JS.jsmodel.x_js[-1] - np.log(args.tract_length), - np.log(args.tract_length) ],
##                           test_JS.x[len(test_JS.jsmodel.x_js):]))
##    test.unpack_x(x)
##
##
##    #test.get_mle()
##    #test.get_individual_summary(summary_file)
##    x_plot = deepcopy(test.x)
##    if not os.path.isdir('./plot'):
##        os.mkdir('./plot')
##
##    if not os.path.isdir('./plot/' + '_'.join(paralog)):
##        os.mkdir('./plot/' + '_'.join(paralog))
##        
##    if args.zoom_in:
##        off_ratio = 0.1
##    else:
##        off_ratio = 0.5
##        
##    tract_length_start = np.exp(-x_js[-1]) * (1.0 - off_ratio)
##    tract_length_end   = np.exp(-x_js[-1]) * (1.0 + off_ratio)
##    tract_length_space = (tract_length_end - tract_length_start)/5
##    tract_length_list = np.log(np.arange(tract_length_start, tract_length_end + tract_length_space, tract_length_space))
##
##    plot_file_name = './plot/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_PSJS_HKY_lnL_'
##    if args.rate_variation:
##        if args.allow_same_codon:
##            plot_file_name += 'rv_SCOK_nonclock_'
##        else:
##            plot_file_name += 'rv_NOSC_nonclock_'
##    else:
##        plot_file_name += 'nonclock_'
##
##    plot_file_name += 'dim_' + str(args.dim) + '_'
##            
##    if args.zoom_in:
##        plot_file_name += 'offratio_' + str(off_ratio) + '_zoomed_in.txt'
##    else:
##        plot_file_name += 'offratio_' + str(off_ratio) +'.txt'
##    with open(plot_file_name, 'w+') as f:
##        tract_length = np.exp(-x_js[-1])
##        #lnL = test.objective_wo_gradient(True, test.x)
##        f.write('\t'.join(['#Tract_length', 'Init_rate', 'lnL', '\n']))
##
##
##
##    for tract_length in tract_length_list:
##        init_rate_list = lnTau - tract_length + np.log(np.arange((1.0 - off_ratio), (1.0 + off_ratio), 2 * off_ratio / 4.0))
##        for init_rate in init_rate_list:
##            x_plot[5] = -tract_length
##            x_plot[4] = init_rate
##            test.unpack_x(x_plot)
##            lnL = test.objective_wo_gradient(True, x_plot)
##            with open(plot_file_name, 'a') as f:
##                f.write('\t'.join([str(np.exp(tract_length)), str(np.exp(init_rate)), str(lnL), '\n']))
