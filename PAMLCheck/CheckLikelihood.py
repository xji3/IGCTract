from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.HMMJSGeneconv import HMMJSGeneconv
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import argparse, os
import numpy as np



if __name__ == '__main__':


    paralog = ['YDR418W', 'YEL054C']

    rate_variation = True

    
    gene_to_orlg_file = '../GeneToOrlg/YDR418W_YEL054C_GeneToOrlg.txt'
    alignment_file = '../MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_input.fasta'
    newicktree = '../YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    Force = None
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    seq_index_file = '../MafftAlignment/YDR418W_YEL054C/YDR418W_YEL054C_seq_index.txt'


###### Now get HKY+PSJS-IGC estimates
    IGC_pm = 'One rate'
    pm_model = 'HKY'
    if not os.path.isdir('./log'):
        os.mkdir('./log')


    if rate_variation:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    else:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])

    

    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                         rate_variation, node_to_pos, terminal_node_list, save_file)
    test_JS.get_mle()

    guess_tract = 10.0
    if rate_variation:
        log_file = './log/PSJS_HKY_rv_guess_' + str(guess_tract) + '_log.txt'
        summary_file = './summary/PSJS_HKY_rv_guess_' + str(guess_tract) + '_summary.txt'
        save_file = './save/PSJS_HKY_rv_guess_' + str(guess_tract) + '_save.txt'
        plot_file = './plot/PSJS_HKY_rv_guess_' + str(guess_tract) + '_plot.txt'
        force = {6:0.0}
    else:
        log_file = './log/PSJS_HKY_guess_' + str(guess_tract) + '_log.txt'
        summary_file = './summary/PSJS_HKY_guess_' + str(guess_tract) + '_summary.txt'
        save_file = './save/PSJS_HKY_guess_' + str(guess_tract) + '_save.txt'
        plot_file = './plot/PSJS_HKY_guess_' + str(guess_tract) + '_plot.txt'
        force = {4:0.0}
    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ]))
    
    PSJS_IGC = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, True, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                      rate_variation, node_to_pos, terminal_node_list, save_file, log_file, force = force)

    x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                           [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ],
                           test_JS.x[len(test_JS.jsmodel.x_js):]))


    x = np.log([0.30297 + 0.22353, 0.30297/(0.30297 + 0.22353), 0.20875/(0.26475+0.20875), 14.463842315785955,\
                0.639751, 11.404187, 1.0, 1.0, \
                0.019194, 0.000004, 0.013513,0.020224,0.001914,0.002955,0.002598,0.006367,0.002979,0.005962,0.005398,0.002868])
    #%AG, %A, %C, kappa, r2, r3
    baseml_edge_to_blen = {
        ('D1', 'N1'):0.000004,
        ('N0', 'D1'):0.0019194,
        ('N0', 'kluyveri'):0.0172746,
        ('N1', 'N2'):0.013513,
        ('N1', 'castellii'):0.020224,
        ('N2', 'N3'):0.001914,
        ('N2', 'bayanus'):0.002955,
        ('N3', 'N4'):0.002598,
        ('N3', 'kudriavzevii'):0.006367,
        ('N4', 'N5'):0.002979,
        ('N4', 'mikatae'):0.005962,
        ('N5', 'cerevisiae'): 0.005398,
        ('N5', 'paradoxus'):0.002868
        }
    PSJS_IGC.unpack_x(x)

    print PSJS_IGC.objective_wo_gradient(True, x)
    

