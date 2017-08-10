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
    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'
    averaged_x = np.loadtxt('./averaged_JS_HKY_One_rate_nonclock_save.txt')
    averaged_x_js, averaged_x_rates = averaged_x[:5], averaged_x[5:]

    save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_' + str(args.guess) + '_nonclock_save.txt'
    log_file  = './log/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_' + str(args.guess) + '_nonclock_log.txt'
    summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_' + str(args.guess) + '_nonclock_summary.txt'

    initial_tract_length_list = np.log([30.0, 200.0])
    guess_lnp = -initial_tract_length_list[args.guess-1]
    
    if args.rate_variation:
        x_js = np.concatenate((averaged_x_js[:-1], np.log([0.7, 3.0]), [averaged_x_js[-1] + guess_lnp, guess_lnp]))
        if args.allow_same_codon:
            save_file = save_file.replace('_nonclock', '_rv_SCOK_nonclock')
            log_file = log_file.replace('_nonclock', '_rv_SCOK_nonclock')
            summary_file = summary_file.replace('_nonclock', '_rv_SCOK_nonclock')
        else:
            save_file = save_file.replace('_nonclock', '_rv_NOSC_nonclock')
            log_file = log_file.replace('_nonclock', '_rv_SCOK_nonclock')
            summary_file = summary_file.replace('_nonclock', '_rv_NOSC_nonclock')
    else:
        x_js = np.concatenate((averaged_x_js[:-1], [averaged_x_js[-1] + guess_lnp, guess_lnp]))

        

     
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, args.cdna, args.allow_same_codon, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                      args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
    #x = np.concatenate((x_js, averaged_x_rates))
    #test.unpack_x(x)


    test.get_mle()
    test.get_individual_summary(summary_file)



if __name__ == '__main__':
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
##    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
##    parser.add_argument('--G', type = int, dest = 'guess', default = 1, help = 'Guess case')
##    parser.add_argument('--heterogeneity', dest = 'rate_variation', action = 'store_true', help = 'rate heterogeneity control')
##    parser.add_argument('--homogeneity', dest = 'rate_variation', action = 'store_false', help = 'rate heterogeneity control')
##    parser.add_argument('--coding', dest = 'cdna', action = 'store_true', help = 'coding sequence control')
##    parser.add_argument('--noncoding', dest = 'cdna', action = 'store_false', help = 'coding sequence control')
##    parser.add_argument('--samecodon', dest = 'allow_same_codon', action = 'store_true', help = 'whether allow pair sites from same codon')
##    parser.add_argument('--no-samecodon', dest = 'allow_same_codon', action = 'store_false', help = 'whether allow pair sites from same codon')
##    
##    main(parser.parse_args())

  

    MyStruct = namedtuple('MyStruct', 'paralog1 paralog2 tract_length dim cdna rate_variation allow_same_codon guess')
    args = MyStruct(paralog1 = 'YLR406C', paralog2 = 'YDL075W', tract_length = 30.0, dim = 1,
                    cdna = True, rate_variation = True, allow_same_codon = True,
                    guess = 1)

    paralog = [args.paralog1, args.paralog2]
    
    gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt'
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_input.fasta'

    tree_newick = '../YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    
    IGC_pm = 'One rate'
    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_seq_index.txt'
    averaged_x = np.loadtxt('./averaged_JS_HKY_One_rate_nonclock_save.txt')
    averaged_x_js, averaged_x_rates = averaged_x[:5], averaged_x[5:]

    save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_' + str(args.guess) + '_nonclock_save.txt'
    log_file  = './log/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_' + str(args.guess) + '_nonclock_log.txt'
    summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_' + str(args.guess) + '_nonclock_summary.txt'

    initial_tract_length_list = np.log([30.0, 200.0])
    guess_lnp = -initial_tract_length_list[args.guess-1]
    
    if args.rate_variation:
        x_js = np.concatenate((averaged_x_js[:-1], np.log([0.7, 3.0]), [averaged_x_js[-1] + guess_lnp, guess_lnp]))
        if args.allow_same_codon:
            save_file = save_file.replace('_nonclock', '_rv_SCOK_nonclock')
            log_file = log_file.replace('_nonclock', '_rv_SCOK_nonclock')
            summary_file = summary_file.replace('_nonclock', '_rv_SCOK_nonclock')
        else:
            save_file = save_file.replace('_nonclock', '_rv_NOSC_nonclock')
            log_file = log_file.replace('_nonclock', '_rv_SCOK_nonclock')
            summary_file = summary_file.replace('_nonclock', '_rv_NOSC_nonclock')
    else:
        x_js = np.concatenate((averaged_x_js[:-1], [averaged_x_js[-1] + guess_lnp, guess_lnp]))

        

     
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, args.cdna, args.allow_same_codon, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                      args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
    #x = np.concatenate((x_js, averaged_x_rates))
    #test.unpack_x(x)


    test.get_mle()
    test.get_individual_summary(summary_file)
