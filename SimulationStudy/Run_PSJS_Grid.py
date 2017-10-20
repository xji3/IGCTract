from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.HMMJSGeneconv import HMMJSGeneconv
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import argparse, os
import numpy as np

def main(args):
    paralog = ['YDR418W', 'YEL054C']
    sim_num = args.sim_num
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
    
    gene_to_orlg_file = '../GeneToOrlg/YDR418W_YEL054C_GeneToOrlg.txt'
    alignment_file = '.' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '.fasta'
    newicktree = './YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    Force = None
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_seq_index.txt'


###### Now get HKY+PSJS-IGC estimates
    IGC_pm = 'One rate'
    pm_model = 'HKY'
    guess_lnP = -np.log(100.0)
    alignment_file = '.' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) + '/YDR418W_YEL054C_sim_' + str(sim_num) + '_newformat.fasta'
    if not os.path.isdir('./log' + case + '' + str(geo) + '_HKY'):
        os.mkdir('./log' + case + '' + str(geo) + '_HKY')
    if not os.path.isdir('./log' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num)):
        os.mkdir('./log' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num))

    if not os.path.isdir('./save' + case + '' + str(geo) + '_HKY'):
        os.mkdir('./save' + case + '' + str(geo) + '_HKY')
    if not os.path.isdir('./save' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num)):
        os.mkdir('./save' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num))

    if not os.path.isdir('./summary' + case + '' + str(geo) + '_HKY'):
        os.mkdir('./summary' + case + '' + str(geo) + '_HKY')
    if not os.path.isdir('./summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num)):
        os.mkdir('./summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num))

    if not os.path.isdir('./plot' + case + '' + str(geo) + '_HKY'):
        os.mkdir('./plot' + case + '' + str(geo) + '_HKY')
    if not os.path.isdir('./plot' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num)):
        os.mkdir('./plot' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num))


    if rate_variation:
        save_file = './save' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = './log' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = './summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    else:
        save_file = './save' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = './log' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = './summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])

    

    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                         rate_variation, node_to_pos, terminal_node_list, save_file)
    test_JS.get_mle()

    
    guess_tract = 1.0
    tract_length_grid = range(5, 65, 5) + range(70, 110, 10) + range(200, 600, 100)
    log_tract_p_list = -np.log(tract_length_grid)
    
    if rate_variation:
        log_file = './log' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' + \
                   str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_nt_log.txt'
        summary_file = './summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' \
                       + str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_nt_summary.txt'
        save_file = './save' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_'\
                    + str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_nt_save.txt'
        plot_file = './plot' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' \
                    + str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_lnL_nt_1D_surface.txt'
        grid_file = './summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_rv_sim_' \
                       + str(sim_num) + '_Tract_' + str(geo) + '_grid_lnL.txt'
    else:
        log_file = './log' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_sim_' \
                   + str(sim_num) + '_Tract_' + str(geo)  + '_guess_' + str(guess_tract) + '_nt_log.txt'
        summary_file = './summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_sim_' \
                       + str(sim_num) + '_Tract_' + str(geo)  + '_guess_' + str(guess_tract) +  '_nt_summary.txt'
        save_file = './save' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_sim_' \
                    + str(sim_num) + '_Tract_' + str(geo)  + '_guess_' + str(guess_tract) +  '_nt_save.txt'
        plot_file = './plot' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_sim_' \
                    + str(sim_num) + '_Tract_' + str(geo)  + '_guess_' + str(guess_tract) +  '_lnL_nt_1D_surface.txt'
        grid_file = './summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/PSJS_HKY_sim_' \
                       + str(sim_num) + '_Tract_' + str(geo) + '_grid_lnL.txt'

    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ]))
    
    PSJS_IGC = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, True, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                      rate_variation, node_to_pos, terminal_node_list, save_file, log_file)

    x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ],
                           test_JS.x[len(test_JS.jsmodel.x_js):]))
    PSJS_IGC.unpack_x(x)
    PSJS_IGC.output_lnL_per_distance_tract_p_list(log_tract_p_list, grid_file)

    


#    log_p_list = np.log(1.0/np.array(range(1, 1001, 2)))
#    PSJS_IGC.plot_tract_p(log_p_list, plot_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--geo', required = True, help = 'Mean tract length')
    parser.add_argument('--sim_num', required = True, help = 'simulation number')
    parser.add_argument('--heterogeneity', dest = 'rate_variation', action = 'store_true', help = 'rate heterogeneity control')
    parser.add_argument('--homogeneity', dest = 'rate_variation', action = 'store_false', help = 'rate heterogeneity control')
    parser.add_argument('--Case', dest = 'Tau_case', default = 'One', help = 'Tau value case')

    
    main(parser.parse_args())


