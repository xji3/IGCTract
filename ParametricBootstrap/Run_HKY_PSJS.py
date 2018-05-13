from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.HMMJSGeneconv import HMMJSGeneconv
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import argparse, os
import numpy as np

def main(args):
    paralog = [args.paralog1, args.paralog2]
    sim_num = args.sim_num
    rate_variation = args.rate_variation
    case_folder = args.case
    
    gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) + '_GeneToOrlg.txt'
    alignment_file = './bootstrapReplicates/' + case_folder + '/' + '_'.join(paralog) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_sim_' + str(sim_num) + '.fasta'
    newicktree = '../YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    Force = None
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_seq_index.txt'


###### Now get HKY+PSJS-IGC estimates
    IGC_pm = 'One rate'
    pm_model = 'HKY'
    log_paralog_folder = './log/' + case_folder + '/' + '_'.join(paralog)
    log_sim_folder = log_paralog_folder + '/sim_' + str(sim_num)
    save_paralog_folder = './save/' + case_folder + '/' + '_'.join(paralog)
    save_sim_folder = save_paralog_folder + '/sim_' + str(sim_num)
    summary_paralog_folder = './summary/' + case_folder + '/' + '_'.join(paralog)
    summary_sim_folder = summary_paralog_folder + '/sim_' + str(sim_num)

    if not os.path.isdir(log_paralog_folder):
        os.mkdir(log_paralog_folder)
    if not os.path.isdir(log_sim_folder):
        os.mkdir(log_sim_folder)
    if not os.path.isdir(save_paralog_folder):
        os.mkdir(save_paralog_folder)
    if not os.path.isdir(save_sim_folder):
        os.mkdir(save_sim_folder)
    if not os.path.isdir(summary_paralog_folder):
        os.mkdir(summary_paralog_folder)
    if not os.path.isdir(summary_sim_folder):
        os.mkdir(summary_sim_folder)


    if rate_variation:
        save_file = save_sim_folder + '/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = log_sim_folder + '/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = summary_sim_folder + '/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    else:
        save_file = save_sim_folder + '/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = log_sim_folder + '/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = summary_sim_folder + '/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])

    

    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                         rate_variation, node_to_pos, terminal_node_list, save_file)
    test_JS.get_mle()
    test_JS.get_individual_summary(summary_file)


    guess_tract = 30.0
    if rate_variation:
        log_file = log_sim_folder + '/PSJS_HKY_rv_sim_' + str(sim_num)  + '_guess_' + str(guess_tract) + '_nt_log.txt'
        summary_file = summary_sim_folder + '/PSJS_HKY_rv_sim_' + str(sim_num) + '_guess_' + str(guess_tract) + '_nt_summary.txt'
        save_file = save_sim_folder + '/PSJS_HKY_rv_sim_' + str(sim_num) + '_guess_' + str(guess_tract) + '_nt_save.txt'
    else:
        log_file = log_sim_folder + '/PSJS_HKY_sim_' + str(sim_num)  + '_guess_' + str(guess_tract) + '_nt_log.txt'
        summary_file = summary_sim_folder + '/PSJS_HKY_sim_' + str(sim_num) + '_guess_' + str(guess_tract) + '_nt_summary.txt'
        save_file = save_sim_folder + '/PSJS_HKY_sim_' + str(sim_num) + '_guess_' + str(guess_tract) + '_nt_save.txt'


    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ]))
    
    PSJS_IGC = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, True, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                      rate_variation, node_to_pos, terminal_node_list, save_file, log_file)

    x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ],
                           test_JS.x[len(test_JS.jsmodel.x_js):]))

    x = np.loadtxt()
    PSJS_IGC.unpack_x(x)

    
    PSJS_IGC.optimize_x_IGC()
    PSJS_IGC.get_individual_summary(summary_file)

#    log_p_list = np.log(1.0/np.array(range(1, 1001, 2)))
#    PSJS_IGC.plot_tract_p(log_p_list, plot_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--sim_num', required = True, help = 'simulation number')
    parser.add_argument('--heterogeneity', dest = 'rate_variation', action = 'store_true', help = 'rate heterogeneity control')
    parser.add_argument('--homogeneity', dest = 'rate_variation', action = 'store_false', help = 'rate heterogeneity control')
    parser.add_argument('--case', dest = 'case', required = True, help = 'which simulated data')

    
    main(parser.parse_args())


