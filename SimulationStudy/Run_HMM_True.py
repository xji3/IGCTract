from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.HMMJSGeneconv import HMMJSGeneconv
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import argparse, os
import numpy as np

def main(args):
    paralog = ['YDR418W', 'YEL054C']
    geo = args.geo
    rate_variation = args.rate_variation
    if args.Tau_case == 'Tenth':
        case = '/TenthTau/Tract_'
    elif args.Tau_case == 'Half':
        case = '/HalfTau/Tract_'
    elif args.Tau_case == 'One':
        case = '/Tract_'
    else:
        raise Exception('Check Tau_case input!')
    
    model = 'HKY'

    for sim_num in range(1, 101):
        alignment_file = '.' + case + '' + str(geo) + '_' + model +'/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '.fasta'
        newicktree = './YeastTree.newick'
        Force = None
        seq_index_file = '../MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_seq_index.txt'

        save_path = './save' + case + '' + str(geo) + '_' + model +'/sim_' + str(sim_num)
        if not os.path.isdir('./save' + case + '' + str(geo)+'_' + model):
            os.mkdir('./save' + case + '' + str(geo)+'_' + model)
        if not os.path.isdir('./save' + case + '' + str(geo) + '_' + model +'/sim_' + str(sim_num)):
            os.mkdir('./save' + case + '' + str(geo) + '_' + model +'/sim_' + str(sim_num))

        if rate_variation:
            save_name = './save' + case + '' + str(geo) + '_' + model +'/sim_' + str(sim_num) + '/YDR418W_YEL054C_Ind_' + model +'_rv_IGC_sim_' + str(sim_num) + '_save.txt'
        else:
            save_name = './save' + case + '' + str(geo) + '_' + model +'/sim_' + str(sim_num) + '/YDR418W_YEL054C_Ind_' + model +'_IGC_sim_' + str(sim_num) + '_save.txt'

        
        Ind_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = None, save_name = save_name,
                                 rate_variation = rate_variation)
        #Ind_IGC.get_mle(True, True, 0, 'BFGS')
        Ind_IGC.update_by_x(np.array([\
        -7.416370120456413639e-01,
        -5.069020547241798180e-01,
        -8.380832454357436401e-01,
        2.671647231636234299e+00,
        -6.176837919230384610e-01,
        2.449282135439945218e+00,
        3.117139303572516784e+00,
        -4.171069746133680312e+00,
        -5.673429216655531349e+00,
        -4.141217212794840563e+00,
        -4.344553499706317545e+00,
        -6.497071863997417651e+00,
        -6.063809695534555289e+00,
        -6.043381749902811961e+00,
        -5.112066696455660697e+00,
        -6.405032887415200271e+00,
        -5.468278611939075162e+00,
        -5.461210564889082519e+00,
        -6.459218907655349895e+00]))

    ###### Now get Ind MG94+IS-IGC+HMM estimates
        summary_path = './summary' + case + '' + str(geo) + '_' + model +'/sim_' + str(sim_num) +'/'
        if not os.path.isdir('./summary' + case + '' + str(geo)+'_' + model):
            os.mkdir('./summary' + case + '' + str(geo)+'_' + model)
        if not os.path.isdir(summary_path):
            os.mkdir(summary_path)
       
        x = np.concatenate((Ind_IGC.x, [0.0]))

        if rate_variation:
            save_file = save_path +  '/True_HMMJS_' + '_'.join(paralog) + '_'+model+'_rv_nonclock_sim_' + str(sim_num) + '_save.txt'
            IGC_sitewise_lnL_file = summary_path + 'True_' +'_'.join(paralog) + '_'+model+'_rv_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'
            NOIGC_sitewise_lnL_file = summary_path + 'True_' + 'NOIGC_' + '_'.join(paralog) + '_'+model+'_rv_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'
        else:
            save_file = save_path +  '/True_HMMJS_' + '_'.join(paralog) + '_'+model+'_nonclock_sim_' + str(sim_num) + '_save.txt'
            IGC_sitewise_lnL_file = summary_path + 'True_' + '_'.join(paralog) + '_'+model+'_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'
            NOIGC_sitewise_lnL_file = summary_path + 'True_' + 'NOIGC_' + '_'.join(paralog) + '_'+model+'_nonclock_sim_' + str(sim_num) + '_sw_lnL.txt'

            
        Ind_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file, False)
        Ind_IGC.get_sitewise_loglikelihood_summary(NOIGC_sitewise_lnL_file, True)

        state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']
         
        HMM_IGC = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,
                         state_list, seq_index_file, model = model, rate_variation = rate_variation)

        if rate_variation:
            summary_file_1D = summary_path + 'True_HMM_' + '_'.join(paralog) + '_' + model + '_rv_nonclock_sim_' + str(sim_num) + '_1D_summary.txt'
            plot_file = './plot' + case + '' + str(geo) + '_' + model + '/sim_' + str(sim_num) + '/True_HMM_' + '_'.join(paralog)+ '_' + model + '_rv_lnL_sim_' + str(sim_num) + '_1D_surface.txt'
        else:
            summary_file_1D = summary_path + 'True_HMM_' + '_'.join(paralog) + '_' + model + '_nonclock_sim_' + str(sim_num) + '_1D_summary.txt'
            plot_file = './plot' + case + '' + str(geo) + '_' + model + '/sim_' + str(sim_num) + '/True_HMM_' + '_'.join(paralog)+ '_' + model + '_lnL_sim_' + str(sim_num) + '_1D_surface.txt'

        HMM_IGC.MG94_IGC.update_by_x(Ind_IGC.x)
        # Plot lnL surface in IS-IGC+HMM model
        log_p_list = np.log(3.0/np.array(range(3, 1001)))
        
        if not os.path.isdir('./plot' + case + '' + str(geo) + '_' + model):
            os.mkdir('./plot' + case + '' + str(geo) + '_' + model)
        if not os.path.isdir('./plot' + case + '' + str(geo) + '_' + model + '/sim_' + str(sim_num)):
            os.mkdir('./plot' + case + '' + str(geo) + '_' + model + '/sim_' + str(sim_num))
        HMM_IGC.plot_tract_p(log_p_list, plot_file)

        # Get MLE and generate summary file
        HMM_IGC.get_mle(display = True, two_step = False, One_Dimension = True)
        HMM_IGC.get_summary(summary_file_1D)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--geo', required = True, help = 'Mean tract length')
    parser.add_argument('--heterogeneity', dest = 'rate_variation', action = 'store_true', help = 'rate heterogeneity control')
    parser.add_argument('--homogeneity', dest = 'rate_variation', action = 'store_false', help = 'rate heterogeneity control')
    parser.add_argument('--Case', dest = 'Tau_case', default = 'One', help = 'Tau value case')
    
    main(parser.parse_args())


##    paralog = ['YDR418W', 'YEL054C']
##    #sim_num = args.sim_num
##    #geo = args.geo
##    #rate_variation = args.rate_variation
##    model = 'HKY'
##    sim_num = 1
##    geo = 3.0
##    rate_variation = True






    