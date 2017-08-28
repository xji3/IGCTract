from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.HMMJSGeneconv import HMMJSGeneconv
import argparse
import numpy as np

def main(args):
    paralog = [args.paralog1, args.paralog2]
    
    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
    newicktree = './YeastTree.newick'
    Force = None
    summary_path = './summary/'
    alignment_file = '../HMMTract/MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    
    Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
    #test.get_mle(True, True, 0, 'BFGS')
    Ind_MG94_IGC.get_sitewise_loglikelihood_summary('./summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt', False)
    Ind_MG94_IGC.get_sitewise_loglikelihood_summary('./summary/NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt', True)
    
    x = np.concatenate((Ind_MG94_IGC.x, [0.0]))
    save_file = './save/HMMJS_' + '_'.join(paralog) + '_MG94_nonclock_save.txt'
    IGC_sitewise_lnL_file = './summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
    NOIGC_sitewise_lnL_file = './summary/NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
    save_path = './save/'
    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_seq_index.txt'
    test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
                     state_list, seq_index_file)
    
    summary_file_1D = './summary/HMM_' + '_'.join(paralog) + '_MG94_nonclock_1D_summary.txt'
    summary_file_all_Dimension = './summary/HMM_' + '_'.join(paralog) + '_MG94_nonclock_all_summary.txt'

    log_p_list = np.log(1.0/np.array(range(1, 1001)))
    plot_file = './plot/HMM_' + '_'.join(paralog) + '_lnL_1D_surface.txt'
    test.plot_tract_p(log_p_list, plot_file)

    
    test.get_mle(display = True, two_step = True, One_Dimension = True)
    test.get_summary(summary_file_1D)
    test.get_mle(display = True, two_step = False, One_Dimension = False)
    test.get_summary(summary_file_all_Dimension)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')

    main(parser.parse_args())

##    pairs = []
##    all_pairs = './Filtered_pairs.txt'
##    with open(all_pairs, 'r') as f:
##        for line in f.readlines():
##            pairs.append(line.replace('\n','').split('_'))
##
##    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
##    newicktree = './YeastTree.newick'
##    Force = None
##    summary_path = './summary/'
##    for pair in pairs[:]:
##        paralog = pair
##        alignment_file = '../HMMTract/MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
##        
##        Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
##        #test.get_mle(True, True, 0, 'BFGS')
##        Ind_MG94_IGC.get_sitewise_loglikelihood_summary('./summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt', False)
##        Ind_MG94_IGC.get_sitewise_loglikelihood_summary('./summary/NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt', True)
##        
##        x = np.concatenate((Ind_MG94_IGC.x, [0.0]))
##        save_file = './save/HMMJS_' + '_'.join(paralog) + '_MG94_nonclock_save.txt'
##        IGC_sitewise_lnL_file = './summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
##        NOIGC_sitewise_lnL_file = './summary/NOIGC_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt'
##        save_path = './save/'
##        seq_index_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_seq_index.txt'
##        test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
##                         state_list, seq_index_file)
##        
##        summary_file_1D = './summary/HMM_' + '_'.join(paralog) + '_MG94_nonclock_1D_summary.txt'
##        summary_file_all_Dimension = './summary/HMM_' + '_'.join(paralog) + '_MG94_nonclock_all_summary.txt'
##        
##        test.get_mle(display = True, two_step = True, One_Dimension = True)
##        test.get_summary(summary_file_1D)
##        test.get_mle(display = True, two_step = False, One_Dimension = False)
##        test.get_summary(summary_file_all_Dimension)
