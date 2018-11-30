from IGCexpansion.CodonSimulator import CodonSimulator
from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
from IGCexpansion.Simulator import Simulator
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
import argparse, os
import numpy as np

if __name__ == '__main__':
    seq_file_list = np.loadtxt('../Harpak/missing_0_species_list.txt', dtype = str)

    pairs = [file_name.split(".")[:2] for file_name in seq_file_list]

    bootstrapReplicateFolder = './bootstrapReplicates'
    if not os.path.isdir(bootstrapReplicateFolder):
        os.mkdir(bootstrapReplicateFolder)

    JSFolder = bootstrapReplicateFolder + '/JSAnalyses'
    if not os.path.isdir(JSFolder):
        os.mkdir(JSFolder)
    PSJSFolder = bootstrapReplicateFolder + '/PSJSAnalyses'
    if not os.path.isdir(PSJSFolder):
        os.mkdir(PSJSFolder)

    gene_to_orlg_file = '../Harpak/GeneToOrlg.txt'
    
    tree_newick = '../Harpak/HarpakTree.newick'
    DupLosList = '../Harpak/HarpakDupLost.txt'
    terminal_node_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    
    IGC_pm = 'One rate'
    cdna = False
    rate_variation = False



    # Simulate according to HKY+rv Model
    display = False
    sim_num = 1
    seed_number_start = 27606

    x_JS = np.loadtxt('./trueValues/JS_HKY_Harpak_One_rate_nonclock_save.txt')
    x_PSJS = np.loadtxt('./trueValues/PSJS_HKY_Harpak_One_rate_nonclock_save.txt')

    for paralog in pairs[-4:]:
        seq_index_file = '../Harpak/prepared_input/' + '_'.join(paralog) +'_seq_index.txt'

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
            x_pm = x_JS[:4]
            rate_variation = False

            x_IGC = [np.exp(x_JS[4]), 1.0]
            init_pm = 'One rate'
            tract_pm = 'One rate'
            pm_IGC = [init_pm, tract_pm]
 
            x_rates = x_JS[5:]
            
            simulator_JS = Simulator(pm_model_name, x_pm, rate_variation,
                             x_IGC, pm_IGC, tree_newick, DupLosList, x_rates,
                             terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file, 
                             cdna)

            print simulator_JS
            simulator_JS.sim(display = display)
            simulator_JS.output_seq(new_format = True)


            save_file = './save/PSJS_HKY_' + '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_1_rv_SCOK_nonclock_save.txt'
            log_file  = './PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_Guess_1_rv_SCOK_nonclock_log.txt'

            x_js = x_PSJS[:6]
            geo = np.exp(-x_PSJS[5])

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
            x_pm = x_PSJS[:4]
            rate_variation = False

            x_IGC = [np.exp(x_PSJS[4]), 1.0 / geo]
            init_pm = 'One rate'
            tract_pm = 'One rate'
            pm_IGC = [init_pm, tract_pm]

            x_rates = x_PSJS[6:]
            
            simulator_PSJS = Simulator(pm_model_name, x_pm, rate_variation,
                             x_IGC, pm_IGC, tree_newick, DupLosList, x_rates,
                             terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file, 
                             cdna)

            print simulator_PSJS
            simulator_PSJS.sim(display = display)
            simulator_PSJS.output_seq(new_format = True)        
