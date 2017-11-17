from IGCexpansion.PSJSGeneconv import *
from IGCexpansion.JSGeneconv import *
import argparse
from collections import namedtuple
import numpy as np

def main(args):
    seq_file = args.seq_file
    
    gene_to_orlg_file = './GeneToOrlg.txt'
    alignment_file = './prepared_input/' + seq_file

    tree_newick = '../HarpakTree.newick'
    DupLosList = '../HarpakDupLost.txt'
    terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    
    IGC_pm = 'One rate'
    seq_index_file = './prepared_input/' + seq_file.replace('.fasta', '_seq_index.txt')

    save_file = './save/' + seq_file.replace('.fasta', '_save.txt')
    log_file = './log/' + seq_file.replace('.fasta', '_log.txt')
    summary_file = './summary/' + seq_file.replace('.fasta', '_summary.txt')

    initial_tract_length_list = np.log([30.0, 200.0, 500.0])
    guess_lnp = -initial_tract_length_list[args.guess-1]
    

    x_js = np.concatenate((np.log([0.3,0.4,0.5,4.0]), [guess_lnp, guess_lnp]))

        

     
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, args.cdna, args.allow_same_codon, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                      args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
    x = np.concatenate((x_js, [np.log(0.1)] * len(averaged_x_rates)))
    test.unpack_x(x)


    test.get_mle(stringent_level = 'high')
    test.get_individual_summary(summary_file)



if __name__ == '__main__':
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--seq', dest= 'seq_file', required = True, help = 'alignment file in fasta format')
##    parser.add_argument('--G', type = int, dest = 'guess', default = 1, help = 'Guess case')
##    
##    main(parser.parse_args())

  

    MyStruct = namedtuple('MyStruct', 'seq_file guess')
    args = MyStruct(seq_file = 'group_317_intron1.fasta',
                    guess = 1)

    seq_file = args.seq_file
    
    gene_to_orlg_file = './GeneToOrlg.txt'
    alignment_file = './prepared_input/' + seq_file

    tree_newick = './HarpakTree.newick'
    DupLosList = './HarpakDupLost.txt'
    terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    
    IGC_pm = 'One rate'
    seq_index_file = './prepared_input/' + seq_file.replace('.fasta', '_seq_index.txt')

    save_file = './save/' + seq_file.replace('.fasta', '_save.txt')
    log_file = './log/' + seq_file.replace('.fasta', '_log.txt')
    summary_file = './summary/' + seq_file.replace('.fasta', '_summary.txt')

    initial_tract_length_list = np.log([30.0, 200.0, 500.0])
    guess_lnp = -initial_tract_length_list[args.guess-1]
    

    x_js = np.concatenate((np.log([0.3,0.4,0.5,4.0]), [guess_lnp, guess_lnp]))

        

     
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, False, False, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                      False, node_to_pos, terminal_node_list, save_file, log_file)



    test.get_mle(stringent_level = 'high')
    test.get_individual_summary(summary_file)
