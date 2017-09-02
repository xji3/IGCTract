from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.HMMJSGeneconv import HMMJSGeneconv
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
import argparse, os
import numpy as np

def main(args):
    paralog = ['YDR418W', 'YEL054C']
    sim_num = args.sim_num
    geo = args.geo
    rate_variation = args.rate_variation

    
    gene_to_orlg_file = '../GeneToOrlg/YDR418W_YEL054C_GeneToOrlg.txt'
    alignment_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '.fasta'
    newicktree = './YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    Force = None
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_seq_index.txt'

    save_path = './save/Tract_' + str(geo) + '/sim_' + str(sim_num)
    if not os.path.isdir('./save/Tract_' + str(geo) ):
        os.mkdir('./save/Tract_' + str(geo))
    if not os.path.isdir('./save/Tract_' + str(geo) + '/sim_' + str(sim_num)):
        os.mkdir('./save/Tract_' + str(geo) + '/sim_' + str(sim_num))
        
    save_name = './save/Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_Ind_MG94_IGC_sim_' + str(sim_num) + '_save.txt'

    
    Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/',
                             save_name = save_name)
    Ind_MG94_IGC.get_mle(True, True, 0, 'BFGS')

###### Now get Ind MG94+IS-IGC+HMM estimates
    summary_path = './summary/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/'
    if not os.path.isdir('./summary/Tract_' + str(geo)):
        os.mkdir('./summary/Tract_' + str(geo))
    if not os.path.isdir(summary_path):
        os.mkdir(summary_path)
   
    x = np.concatenate((Ind_MG94_IGC.x, [0.0]))
    save_file = save_path +  '/HMMJS_' + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_save.txt'
    IGC_sitewise_lnL_file = summary_path + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'
    NOIGC_sitewise_lnL_file = summary_path + 'NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'
        
    Ind_MG94_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file, False)
    Ind_MG94_IGC.get_sitewise_loglikelihood_summary(NOIGC_sitewise_lnL_file, True)

    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
     
    HMM_MG94_IGC = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
                     state_list, seq_index_file)

    summary_file_1D = summary_path + 'HMM_' + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_1D_summary.txt'
    
    # Plot lnL surface in IS-IGC+HMM model
    log_p_list = np.log(3.0/np.array(range(3, 1001)))
    plot_file = './plot/Tract_' + str(geo) + '/sim_' + str(sim_num) + '/HMM_' + '_'.join(paralog) + '_lnL_sim_' + str(sim_num) + '_1D_surface.txt'
    if not os.path.isdir('./plot/Tract_' + str(geo)):
        os.mkdir('./plot/Tract_' + str(geo))
    if not os.path.isdir('./plot/Tract_' + str(geo) + '/sim_' + str(sim_num)):
        os.mkdir('./plot/Tract_' + str(geo) + '/sim_' + str(sim_num))
    HMM_MG94_IGC.plot_tract_p(log_p_list, plot_file)

    # Get MLE and generate summary file
    HMM_MG94_IGC.get_mle(display = True, two_step = True, One_Dimension = True)
    HMM_MG94_IGC.get_summary(summary_file_1D)

###### Now get HKY+PSJS-IGC estimates
    IGC_pm = 'One rate'
    pm_model = 'HKY'
    guess_lnP = -np.log(100.0)
    alignment_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '_newformat.fasta'
    if not os.path.isdir('./log/Tract_' + str(geo)):
        os.mkdir('./log/Tract_' + str(geo))
    if not os.path.isdir('./log/Tract_' + str(geo) + '/sim_' + str(sim_num)):
        os.mkdir('./log/Tract_' + str(geo) + '/sim_' + str(sim_num))
        
    if rate_variation:
        x_js = np.concatenate((Ind_MG94_IGC.x_process[:-2], np.log([0.7, 3.0]), [Ind_MG94_IGC.x_process[-1] + guess_lnP, guess_lnP]))
        log_file = './log/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
        summary_file = './summary/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
        save_file = './save/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
        plot_file = './plot/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_lnL_1D_surface.txt'
    else:
        x_js = np.concatenate((Ind_MG94_IGC.x_process[:-2], [Ind_MG94_IGC.x_process[-1] + guess_lnP, guess_lnP]))
        log_file = './log/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
        summary_file = './summary/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
        save_file = './save/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
        plot_file = './plot/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_lnL_1D_surface.txt'

    PSJS_IGC = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, True, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                      rate_variation, node_to_pos, terminal_node_list, save_file, log_file)

    PSJS_IGC.optimize_x_IGC()
    PSJS_IGC.get_individual_summary(summary_file)
    PSJS_IGC.plot_tract_p(log_p_list, plot_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--geo', required = True, help = 'Mean tract length')
    parser.add_argument('--sim_num', required = True, help = 'simulation number')
    parser.add_argument('--heterogeneity', dest = 'rate_variation', action = 'store_true', help = 'rate heterogeneity control')
    parser.add_argument('--homogeneity', dest = 'rate_variation', action = 'store_false', help = 'rate heterogeneity control')
    
    main(parser.parse_args())


######## Get MLE of IS-IGC model
##    paralog = ['YDR418W', 'YEL054C']
##    #sim_num = args.sim_num
##    #geo = args.geo
##    #rate_variation = args.rate_variation
##    sim_num = 1
##    geo = 3.0
##    rate_variation = True
##    
##    gene_to_orlg_file = '../GeneToOrlg/YDR418W_YEL054C_GeneToOrlg.txt'
##    alignment_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '.fasta'
##    newicktree = './YeastTree.newick'
##    DupLosList = '../YeastTestDupLost.txt'
##    Force = None
##    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
##    node_to_pos = {'D1':0}
##    seq_index_file = '../MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_seq_index.txt'
##
##    save_path = './save/Tract_' + str(geo) + '/sim_' + str(sim_num)
##    if not os.path.isdir('./save/Tract_' + str(geo) ):
##        os.mkdir('./save/Tract_' + str(geo))
##    if not os.path.isdir('./save/Tract_' + str(geo) + '/sim_' + str(sim_num)):
##        os.mkdir('./save/Tract_' + str(geo) + '/sim_' + str(sim_num))
##        
##    save_name = './save/Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_Ind_MG94_IGC_sim_' + str(sim_num) + '_save.txt'
##
##    
##    Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/',
##                             save_name = save_name)
##    Ind_MG94_IGC.get_mle(True, True, 0, 'BFGS')
##
######## Now get Ind MG94+IS-IGC+HMM estimates
##    summary_path = './summary/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/'
##    if not os.path.isdir('./summary/Tract_' + str(geo)):
##        os.mkdir('./summary/Tract_' + str(geo))
##    if not os.path.isdir(summary_path):
##        os.mkdir(summary_path)
##   
##    x = np.concatenate((Ind_MG94_IGC.x, [0.0]))
##    save_file = save_path +  '/HMMJS_' + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_save.txt'
##    IGC_sitewise_lnL_file = summary_path + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'
##    NOIGC_sitewise_lnL_file = summary_path + 'NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'
##        
##    Ind_MG94_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file, False)
##    Ind_MG94_IGC.get_sitewise_loglikelihood_summary(NOIGC_sitewise_lnL_file, True)
##
##    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
##     
##    HMM_MG94_IGC = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
##                     state_list, seq_index_file)
##
##    summary_file_1D = summary_path + 'HMM_' + '_'.join(paralog) + '_MG94_nonclock_sim_' + str(sim_num) + '_1D_summary.txt'
##    
##    # Plot lnL surface in IS-IGC+HMM model
##    log_p_list = np.log(3.0/np.array(range(3, 1001)))
##    plot_file = './plot/Tract_' + str(geo) + '/sim_' + str(sim_num) + '/HMM_' + '_'.join(paralog) + '_lnL_sim_' + str(sim_num) + '_1D_surface.txt'
##    if not os.path.isdir('./plot/Tract_' + str(geo)):
##        os.mkdir('./plot/Tract_' + str(geo))
##    if not os.path.isdir('./plot/Tract_' + str(geo) + '/sim_' + str(sim_num)):
##        os.mkdir('./plot/Tract_' + str(geo) + '/sim_' + str(sim_num))
##    HMM_MG94_IGC.plot_tract_p(log_p_list, plot_file)
##
##    # Get MLE and generate summary file
##    HMM_MG94_IGC.get_mle(display = True, two_step = True, One_Dimension = True)
##    HMM_MG94_IGC.get_summary(summary_file_1D)
##
######## Now get HKY+PSJS-IGC estimates
##    IGC_pm = 'One rate'
##    pm_model = 'HKY'
##    guess_lnP = -np.log(100.0)
##    alignment_file = './Tract_' + str(geo) + '/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '_newformat.fasta'
##    if not os.path.isdir('./log/Tract_' + str(geo)):
##        os.mkdir('./log/Tract_' + str(geo))
##    if not os.path.isdir('./log/Tract_' + str(geo) + '/sim_' + str(sim_num)):
##        os.mkdir('./log/Tract_' + str(geo) + '/sim_' + str(sim_num))
##        
##    if rate_variation:
##        x_js = np.concatenate((Ind_MG94_IGC.x_process[:-2], np.log([0.7, 3.0]), [Ind_MG94_IGC.x_process[-1] + guess_lnP, guess_lnP]))
##        log_file = './log/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
##        summary_file = './summary/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
##        save_file = './save/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
##        plot_file = './plot/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_lnL_1D_surface.txt'
##    else:
##        x_js = np.concatenate((Ind_MG94_IGC.x_process[:-2], [Ind_MG94_IGC.x_process[-1] + guess_lnP, guess_lnP]))
##        log_file = './log/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
##        summary_file = './summary/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
##        save_file = './save/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_log.txt'
##        plot_file = './plot/Tract_' + str(geo) + '/sim_' + str(sim_num) +'/PSJS_HKY_sim_' + str(sim_num) + '_Tract_' + str(int(geo)) + '_lnL_1D_surface.txt'
##
##    PSJS_IGC = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, True, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
##                      rate_variation, node_to_pos, terminal_node_list, save_file, log_file)
##
##    PSJS_IGC.optimize_x_IGC()
##    PSJS_IGC.get_individual_summary(summary_file)
##    PSJS_IGC.plot_tract_p(log_p_list, plot_file)
    
