from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.PSJSGeneconv import PSJSGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
import argparse, os
import numpy as np

def main(args):
    paralog = ['YDR418W', 'YEL054C']
    div = int(args.div)
    for sim_num in range(div*10 + 1, (div+1)*10 + 1):
        geo = 100.0
        rate_variation = True
        case = '/Tract_'
        
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
    ##    test_JS.unpack_x(np.array([\
    ##        -7.416370120456413639e-01,
    ##        -5.069020547241798180e-01,
    ##        -8.380832454357436401e-01,
    ##        2.671647231636234299e+00,
    ##        -6.176837919230384610e-01,
    ##        2.449282135439945218e+00,
    ##        3.117139303572516784e+00,
    ##        -4.171069746133680312e+00,
    ##        -5.673429216655531349e+00,
    ##        -4.141217212794840563e+00,
    ##        -4.344553499706317545e+00,
    ##        -6.497071863997417651e+00,
    ##        -6.063809695534555289e+00,
    ##        -6.043381749902811961e+00,
    ##        -5.112066696455660697e+00,
    ##        -6.405032887415200271e+00,
    ##        -5.468278611939075162e+00,
    ##        -5.461210564889082519e+00,
    ##        -6.459218907655349895e+00]))
        test_JS.get_mle()

        guess_tract = 20.0

        log_file = './log' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/True_PSJS_HKY_rv_sim_' + \
                   str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_nt_log.txt'
        summary_file = './summary' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/True_PSJS_HKY_rv_sim_' \
                       + str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_nt_summary.txt'
        save_file = './save' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/True_PSJS_HKY_rv_sim_'\
                    + str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_nt_save.txt'
        plot_file = './plot' + case + '' + str(geo) + '_HKY/sim_' + str(sim_num) +'/True_PSJS_HKY_rv_sim_' \
                    + str(sim_num) + '_Tract_' + str(geo) + '_guess_' + str(guess_tract) + '_lnL_nt_1D_surface.txt'

        x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                               [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ]))
        
        PSJS_IGC = PSJSGeneconv(alignment_file, gene_to_orlg_file, seq_index_file, True, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                          rate_variation, node_to_pos, terminal_node_list, save_file, log_file)

        x = np.concatenate((test_JS.jsmodel.x_js[:-1], \
                               [ test_JS.jsmodel.x_js[-1] - np.log(guess_tract), - np.log(guess_tract) ],
                               test_JS.x[len(test_JS.jsmodel.x_js):]))
        PSJS_IGC.unpack_x(x)

        
        PSJS_IGC.get_mle()
        PSJS_IGC.get_individual_summary(summary_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--div', dest = 'div', help = '1 of 10 simulation batches')

    main(parser.parse_args())


