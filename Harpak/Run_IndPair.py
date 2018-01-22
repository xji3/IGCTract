# This is a separate file to run all 20 Harpak data sets together
# This file implements a naiive parallel computing method without careful testing
# Xiang Ji
# xji3@ncsu.edu


# OK, naiive parallel computing upgrade
# http://sebastianraschka.com/Articles/2014_multiprocessing.html

from IGCexpansion.PSJSGeneconv import PSJSGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import numpy as np
import scipy, os, argparse
import multiprocessing as mp
from functools import partial
from collections import namedtuple
from Run_JS import Run_JS_Harpak_all, Run_PSJS_Harpak_all


def main(args):
    seq_file_list = np.loadtxt('missing_0_species_list.txt', dtype = str)
    # Now group alignment seq files by gene pair
    group_name_list = sorted(list(set([seq_name[:9] for seq_name in seq_file_list])))

    group_to_seq_file_dict = {group_name:[seq_file for seq_file in seq_file_list if seq_file[:9] == group_name] for group_name in group_name_list}

    group_name = group_name_list[args.group]

    alignment_file_list = ['./prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                          + '.fasta' for seq_file in group_to_seq_file_dict[group_name]]
    seq_index_file_list = ['./prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                          + '_seq_index.txt' for seq_file in group_to_seq_file_dict[group_name]]

    save_file_list = ['./save/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                          + '_guess_' + str(args.guess) + '_IndGroup_save.txt' for seq_file in group_to_seq_file_dict[group_name]]
    log_file_list = ['./log/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                          + '_guess_' + str(args.guess) + '_IndGroup_log.txt' for seq_file in group_to_seq_file_dict[group_name]]
    summary_file_list = ['./summary/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                          + '_guess_' + str(args.guess) + '_IndGroup_summary.txt' for seq_file in group_to_seq_file_dict[group_name]]
    gradient_file_list = ['./summary/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                          + '_guess_' + str(args.guess) + '_IndGroup_gradient.txt' for seq_file in group_to_seq_file_dict[group_name]]
    hessian_file_list = ['./summary/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
                          + '_guess_' + str(args.guess) + '_IndGroup_essian.txt' for seq_file in group_to_seq_file_dict[group_name]]

    gene_to_orlg_file = './GeneToOrlg.txt'
    save_file = './save/PSJS_' +  group_name + '_Grand_save_guess_' + str(args.guess) + '.txt'

    tree_newick = './HarpakTree.newick'
    DupLosList = './HarpakDupLost.txt'
    terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'

    IGC_pm = 'One rate'

    initial_tract_length_list = np.log([30.0, 200.0, 500.0])
    guess_lnp = -initial_tract_length_list[args.guess - 1]

    rate_variation = False
    x_js = np.concatenate((np.log([0.3,0.4,0.5,4.0]), [guess_lnp, guess_lnp]))


    # multiprocess_combined_list = [[0, 9, 12], [7, 5], [6, 1], [18, 8],
    #                              [2], [3], [4], [10], [11],
    #                              [13], [14], [15], [16], [17],
    #                              [19]]
     
    test = Run_PSJS_Harpak_all(alignment_file_list, gene_to_orlg_file,
                              seq_index_file_list,
                              tree_newick, DupLosList,
                              x_js, pm_model, IGC_pm, rate_variation,
                              node_to_pos, terminal_node_list,
                              save_file_list, log_file_list, save_file
                              )
    test.get_mle(stringent_level = 'high')
    Godambe_x = np.array([test.psjsgeneconv_list[0].psjsmodel.x_IGC[0] - test.psjsgeneconv_list[0].psjsmodel.x_IGC[1],\
                         test.psjsgeneconv_list[0].psjsmodel.x_IGC[1] - np.log(1.0 - np.exp(test.psjsgeneconv_list[0].psjsmodel.x_IGC[1]))])
    test.get_gradient_hessian(Godambe_x, gradient_file_list, hessian_file_list)
   
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--Guess', type = int, dest = 'guess', default = 1, help = 'Guess case')
    parser.add_argument('--Group', type = int, dest = 'group', default = 1, help = 'Group number starting from 0')

    main(parser.parse_args())

    # MyStruct = namedtuple('MyStruct', 'guess group')
    # args = MyStruct(guess = 1, group = 0)

    # seq_file_list = np.loadtxt('missing_0_species_list.txt', dtype = str)
    # # Now group alignment seq files by gene pair
    # group_name_list = sorted(list(set([seq_name[:9] for seq_name in seq_file_list])))

    # group_to_seq_file_dict = {group_name:[seq_file for seq_file in seq_file_list if seq_file[:9] == group_name] for group_name in group_name_list}

    # group_name = group_name_list[args.group]

    # alignment_file_list = ['./prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
    #                       + '.fasta' for seq_file in group_to_seq_file_dict[group_name]]
    # seq_index_file_list = ['./prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
    #                       + '_seq_index.txt' for seq_file in group_to_seq_file_dict[group_name]]

    # save_file_list = ['./save/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
    #                       + '_guess_' + str(args.guess) + '_IndGroup_save.txt' for seq_file in group_to_seq_file_dict[group_name]]
    # log_file_list = ['./log/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
    #                       + '_guess_' + str(args.guess) + '_IndGroup_log.txt' for seq_file in group_to_seq_file_dict[group_name]]
    # summary_file_list = ['./summary/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
    #                       + '_guess_' + str(args.guess) + '_IndGroup_summary.txt' for seq_file in group_to_seq_file_dict[group_name]]
    # gradient_file_list = ['./summary/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
    #                       + '_guess_' + str(args.guess) + '_IndGroup_gradient.txt' for seq_file in group_to_seq_file_dict[group_name]]
    # hessian_file_list = ['./summary/PSJS_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
    #                       + '_guess_' + str(args.guess) + '_IndGroup_essian.txt' for seq_file in group_to_seq_file_dict[group_name]]

    # gene_to_orlg_file = './GeneToOrlg.txt'
    # save_file = './save/PSJS_' +  group_name + '_Grand_save_guess_' + str(args.guess) + '.txt'

    # tree_newick = './HarpakTree.newick'
    # DupLosList = './HarpakDupLost.txt'
    # terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
    # node_to_pos = {'D1':0}
    # pm_model = 'HKY'

    # IGC_pm = 'One rate'

    # initial_tract_length_list = np.log([30.0, 200.0, 500.0])
    # guess_lnp = -initial_tract_length_list[args.guess - 1]

    # rate_variation = False
    # x_js = np.concatenate((np.log([0.3,0.4,0.5,4.0]), [guess_lnp, guess_lnp]))


    # # multiprocess_combined_list = [[0, 9, 12], [7, 5], [6, 1], [18, 8],
    # #                              [2], [3], [4], [10], [11],
    # #                              [13], [14], [15], [16], [17],
    # #                              [19]]
     
    # test = Run_PSJS_Harpak_all(alignment_file_list, gene_to_orlg_file,
    #                           seq_index_file_list,
    #                           tree_newick, DupLosList,
    #                           x_js, pm_model, IGC_pm, rate_variation,
    #                           node_to_pos, terminal_node_list,
    #                           save_file_list, log_file_list, save_file
    #                           )
    # test.get_mle(stringent_level = 'high')
    # Godambe_x = np.array([test.psjsgeneconv_list[0].psjsmodel.x_IGC[0] - test.psjsgeneconv_list[0].psjsmodel.x_IGC[1],\
    #                      test.psjsgeneconv_list[0].psjsmodel.x_IGC[1] - np.log(1.0 - np.exp(test.psjsgeneconv_list[0].psjsmodel.x_IGC[1]))])
    # test.get_gradient_hessian(Godambe_x, gradient_file_list, hessian_file_list)
   
   

# #######################################################################################################################################
# ###############################           IS-IGC results              #################################################################
# #######################################################################################################################################
# ####
#     seq_file_list = np.loadtxt('missing_0_species_list.txt', dtype = str)
#     alignment_file_list = ['./prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                            + '.fasta' for seq_file in seq_file_list]
#     seq_index_file_list = ['./prepared_input/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                            + '_seq_index.txt' for seq_file in seq_file_list]

#     save_file_list = ['./save/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                            + '_guess_' + str(args.guess) + '_HKY_JS_save.txt' for seq_file in seq_file_list]
#     log_file_list = ['./log/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                            + '_guess_' + str(args.guess) + '_HKY_JS_log.txt' for seq_file in seq_file_list]
#     summary_file_list = ['./summary/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                            + '_guess_' + str(args.guess) + '_HKY_JS_summary.txt' for seq_file in seq_file_list]

#     gradient_file_list = ['./summary/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                           + '_guess_' + str(args.guess) + '_HKY_JS_gradient.txt' for seq_file in seq_file_list]

#     hessian_file_list = ['./summary/' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                           + '_guess_' + str(args.guess) + '_HKY_JS_hessian.txt' for seq_file in seq_file_list]

#     gene_to_orlg_file = './GeneToOrlg.txt'
#     save_file = './save/Grand_save_HKY_JS_guess_' + str(args.guess) + '.txt'

#     tree_newick = './HarpakTree.newick'
#     DupLosList = './HarpakDupLost.txt'
#     terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
#     node_to_pos = {'D1':0}
#     pm_model = 'HKY'
    
#     IGC_pm = 'One rate'

#     initial_tract_length_list = np.log([30.0, 200.0, 500.0])
#     guess_lnp = -initial_tract_length_list[args.guess - 1]

#     rate_variation = False
    

#     x_js = np.log([0.3,0.4,0.5,4.0, 3.0])


# ##    alignment_file_list = alignment_file_list[:2]
# ##    seq_index_file_list = seq_index_file_list[:2]
# ##    save_file_list = save_file_list[:2]
# ##    log_file_list = log_file_list[:2]
      
#     test = Run_JS_Harpak_all(alignment_file_list, gene_to_orlg_file,
#                  seq_index_file_list, 
#                  tree_newick, DupLosList,
#                  x_js, pm_model, IGC_pm, rate_variation,
#                  node_to_pos, terminal_node_list,
#                  save_file_list, log_file_list, save_file)

#     self = test

#     #results = test.objective_and_gradient(True, test.x)
#     #print results
#     # test.get_mle(stringent_level = 'high')
#     for js in test.jsgeneconv_list:
#         js.get_expectedNumGeneconv()
#         js.get_expectedMutationNum()
#     test.get_summary(summary_file_list)
#     test.get_Godambe_matrix(test.x, gradient_file_list, hessian_file_list)
    
#     force = {4:0.0}

#     save_file_list = ['./save/Force_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                            + '_guess_' + str(args.guess) + '_HKY_JS_save.txt' for seq_file in seq_file_list]
#     log_file_list = ['./log/Force_' + seq_file.replace('.pos.seq.formatted', '').replace('.', '_') \
#                            + '_guess_' + str(args.guess) + '_HKY_JS_log.txt' for seq_file in seq_file_list]
#     save_file = './save/Force_Grand_save_HKY_JS_guess_' + str(args.guess) + '.txt'

    
#     JS_IGC_Force = Run_JS_Harpak_all(alignment_file_list, gene_to_orlg_file,
#                  seq_index_file_list, 
#                  tree_newick, DupLosList,
#                  x_js, pm_model, IGC_pm, rate_variation,
#                  node_to_pos, terminal_node_list,
#                  save_file_list, log_file_list, save_file, force = force)
#     JS_IGC_Force.get_mle(stringent_level = 'high')


                          
