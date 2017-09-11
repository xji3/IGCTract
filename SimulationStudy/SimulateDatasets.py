from IGCexpansion.CodonSimulator import CodonSimulator
from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
import argparse, os
import numpy as np

if __name__ == '__main__':

###### Check MLE
    paralog = ['YDR418W', 'YEL054C']
    gene_to_orlg_file = '../GeneToOrlg/YDR418W_YEL054C_GeneToOrlg.txt'
    alignment_file = '../HMMTract/MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_input.fasta'
    newicktree = './YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    Force = None
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_seq_index.txt'

    Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/',
                             save_name = './save/YDR418W_YEL054C_Ind_MG94_IGC_save.txt')
    #Ind_MG94_IGC.get_mle(True, True, 0, 'BFGS')

###### Now simulate datasets
    display = False
    mean_tract_list = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]
    sim_num = 1
    seed_number_start = 27606
    
    for geo in mean_tract_list[:]:
        for sim_num in range(1, 101):

            seed_number = seed_number_start + sim_num

            if not os.path.isdir('./Tract_' + str(geo)):
                os.mkdir('./Tract_' + str(geo))
            if not os.path.isdir('./Tract_' + str(geo) + '/sim_' + str(sim_num)):
                os.mkdir('./Tract_' + str(geo) + '/sim_' + str(sim_num))
            seq_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '.fasta'
            IGC_log_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '_IGC.log'
            PM_log_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '_PM.log'

            pm_model_name = 'MG94'
            #x_pm = np.log([0.4, 0.5, 0.2, 9.2, 1.0])
            x_pm = Ind_MG94_IGC.x_process[:-1]
            rate_variation = False

            x_IGC = [Ind_MG94_IGC.tau * 3.0 / geo, 3.0 / geo]
            init_pm = 'One rate'
            tract_pm = 'One rate'
            pm_IGC = [init_pm, tract_pm]

        ##            x_rates = [-4.170654939766711422e+00,
        ##                       -5.674236262981605883e+00,
        ##                       -4.140979602575983520e+00,
        ##                       -4.344239699023852097e+00,
        ##                       -6.496123290482403334e+00,
        ##                       -6.063647134296714647e+00,
        ##                       -6.043806966727234276e+00,
        ##                       -5.111657692573940537e+00,
        ##                       -6.404488905061815451e+00,
        ##                       -5.467996717925044159e+00,
        ##                       -5.460686727891754799e+00,
        ##                       -6.459940982759793116e+00]
            
            x_rates = Ind_MG94_IGC.x_rates
            
            test = CodonSimulator(pm_model_name, x_pm, rate_variation,
                             x_IGC, pm_IGC, newicktree, DupLosList, x_rates,
                             terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file)

            self = test
            print test
            
            test.sim_root()
            #edge = ('N0', 'kluyveri')
            #test.sim_one_branch(edge, True)

            test.sim(display = display)
            test.output_seq()

            test.seq_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '_newformat.fasta'
            test.output_seq(True)
