from IGCexpansion.CodonSimulator import CodonSimulator
from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
from IGCexpansion.Simulator import Simulator
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
import argparse, os
import numpy as np

if __name__ == '__main__':

    pairs = []
    all_pairs = '../Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
       for line in f.readlines():
           pairs.append(line.replace('\n','').split('_'))

    bootstrapReplicateFolder = './bootstrapReplicates'
    if not os.path.isdir(bootstrapReplicateFolder):
        os.mkdir(bootstrapReplicateFolder)

    JSFolder = bootstrapReplicateFolder + '/JSAnalyses'
    if not os.path.isdir(JSFolder):
        os.mkdir(JSFolder)
    PSJSFolder = bootstrapReplicateFolder + '/PSJSAnalyses'
    if not os.path.isdir(PSJSFolder):
        os.mkdir(PSJSFolder)
    

    paralog = ["EDN", "ECP"]

    gene_to_orlg_file = './' + '_'.join(paralog) +'_GeneToOrlg.txt'
    alignment_file = './' + '_'.join(paralog) +'_Cleaned_input.fasta'

    tree_newick = './input_tree.newick'
    DupLosList = './EDN_ECP_DupLost.txt'
    terminal_node_list = ['Tamarin', 'Macaque', 'Orangutan', 'Gorilla', 'Chimpanzee']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    cdna = True
    rate_variation = True
    
    IGC_pm = 'One rate'
    seq_index_file = './' + '_'.join(paralog) +'_seq_index.txt'

    save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') +  '_nonclock_save.txt'
    log_file  = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
    summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
 
    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    #force = {6:0.0}

    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                    rate_variation, node_to_pos, terminal_node_list, save_file)
    test_JS.get_mle()
    test_JS.get_individual_summary(summary_file)

    # Simulate according to HKY+rv Model
    display = False
    sim_num = 1
    seed_number_start = 27606

    for sim_num in range(1,101):
        geo = 1
        seed_number = seed_number_start + sim_num
        paralog_folder = './bootstrapReplicates/JSAnalyses/' + '_'.join(paralog) 
        sim_folder = './bootstrapReplicates/JSAnalyses/' + '_'.join(paralog) + '/sim_' + str(sim_num)
        if not os.path.isdir(paralog_folder):
            os.mkdir(paralog_folder)
        if not os.path.isdir(sim_folder):
            os.mkdir(sim_folder)

        seq_file = sim_folder + '/'  + '_'.join(paralog) + '_sim_' + str(sim_num) + '.fasta'
        IGC_log_file = sim_folder + '/'  + '_'.join(paralog) + '_sim_' + str(sim_num) + '_IGC.log'
        PM_log_file = sim_folder + '/'  + '_'.join(paralog) + '_sim_' + str(sim_num) + '_PM.log'

        pm_model_name = 'HKY'
        x_pm = test_JS.jsmodel.x_pm
        rate_variation = True

        x_IGC = [test_JS.jsmodel.IGCModel.parameters['Tau'] / geo, 1.0 / geo]
        init_pm = 'One rate'
        tract_pm = 'One rate'
        pm_IGC = [init_pm, tract_pm]


        if test_JS.root_by_dup:
            x_rates = test_JS.x[-len(test_JS.tree.edge_list):]
        else:
            x_rates = test_JS.x[-(len(test_JS.tree.edge_list) - 1):]
        
        simulator_JS = Simulator(pm_model_name, x_pm, rate_variation,
                         x_IGC, pm_IGC, tree_newick, DupLosList, x_rates,
                         terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file)

        print simulator_JS
        simulator_JS.sim(display = display)
        simulator_JS.output_seq(new_format = True)


        save_file = './save/PSJS_HKY_' + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_1_rv_SCOK_nonclock_save.txt'
        log_file  = './PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_1_rv_SCOK_nonclock_log.txt'

        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3, 0.4])
        test_PSJS = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, True, True, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                      rate_variation, node_to_pos, terminal_node_list, save_file, log_file)

        geo = np.exp(-test_PSJS.psjsmodel.x_IGC[1])

        paralog_folder = './bootstrapReplicates/PSJSAnalyses/' + '_'.join(paralog) 
        sim_folder = './bootstrapReplicates/PSJSAnalyses/' + '_'.join(paralog) + '/sim_' + str(sim_num)
        if not os.path.isdir(paralog_folder):
            os.mkdir(paralog_folder)
        if not os.path.isdir(sim_folder):
            os.mkdir(sim_folder)

        seq_file = sim_folder + '/'  + '_'.join(paralog) + '_sim_' + str(sim_num) + '.fasta'
        IGC_log_file = sim_folder + '/'  + '_'.join(paralog) + '_sim_' + str(sim_num) + '_IGC.log'
        PM_log_file = sim_folder + '/'  + '_'.join(paralog) + '_sim_' + str(sim_num) + '_PM.log'

        pm_model_name = 'HKY'
        x_pm = test_PSJS.psjsmodel.x_pm
        rate_variation = True

        x_IGC = [np.exp(test_PSJS.psjsmodel.x_IGC[0]), 1.0 / geo]
        init_pm = 'One rate'
        tract_pm = 'One rate'
        pm_IGC = [init_pm, tract_pm]


        if test_JS.root_by_dup:
            x_rates = test_PSJS.x[-len(test_PSJS.tree.edge_list):]
        else:
            x_rates = test_PSJS.x[-(len(test_PSJS.tree.edge_list) - 1):]
        
        simulator_PSJS = Simulator(pm_model_name, x_pm, rate_variation,
                         x_IGC, pm_IGC, tree_newick, DupLosList, x_rates,
                         terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file)

        print simulator_PSJS
        simulator_PSJS.sim(display = display)
        simulator_PSJS.output_seq(new_format = True)        
