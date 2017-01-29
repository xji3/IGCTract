from IGCexpansion.PSJSGeneconv import *
from IGCexpansion.JSGeneconv import *
import argparse
from collections import namedtuple

def main(args):
    paralog = [args.paralog1, args.paralog2]
    
    gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt'
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_input.fasta'

    tree_newick = '../YeastTree.newick'
    DupLosList = '../YeastTestDupLost.txt'
    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])
    IGC_pm = 'One rate'

    save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
    summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'

    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
                      node_to_pos, terminal_node_list, save_file)
    test_JS.get_mle()
    test_JS.get_expectedNumGeneconv()
    test_JS.get_individual_summary(summary_file)
    
    save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
    summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'

    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], [-0.6, -0.6 - test_JS.jsmodel.x_js[-1]] ))
    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
                      node_to_pos, terminal_node_list, save_file)
    test.get_mle()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    
    main(parser.parse_args())
##
##    MyStruct = namedtuple('MyStruct', 'paralog1 paralog2')
##    args = MyStruct(paralog1 = 'YDR418W', paralog2 = 'YEL054C')
##
##    paralog = [args.paralog1, args.paralog2]
##    
##    gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt'
##    alignment_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_input.fasta'
##
##    tree_newick = '../YeastTree.newick'
##    DupLosList = '../YeastTestDupLost.txt'
##    terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
##    node_to_pos = {'D1':0}
##    pm_model = 'HKY'
##    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])
##    IGC_pm = 'One rate'
##
##    save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
##    summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
##
##    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
##                      node_to_pos, terminal_node_list, save_file)
##    test_JS.get_mle()
##    test_JS.get_expectedNumGeneconv()
##    test_JS.get_individual_summary(summary_file)
##    
##    save_file = './save/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
##    summary_file = './summary/PSJS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
##
##    x_js = np.concatenate((test_JS.jsmodel.x_js[:-1], [-0.6, -0.6 - test_JS.jsmodel.x_js[-1]] ))
##    test = PSJSGeneconv(alignment_file, gene_to_orlg_file, tree_newick, DupLosList,x_js, pm_model, IGC_pm,
##                      node_to_pos, terminal_node_list, save_file)
##    test.get_mle()
