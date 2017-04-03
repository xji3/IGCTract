# This script is used to get lnL summary of JS models
# Xiang Ji
# xji3@ncsu.edu

from IGCexpansion.JSGeneconv import *

if __name__ == '__main__':
    pairs = []
    all_pairs = '../Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    cdna = True
    for rate_variation in [True, False]:

        for paralog in pairs:

            gene_to_orlg_file = '../GeneToOrlg/' + '_'.join(paralog) +'_GeneToOrlg.txt'
            alignment_file = '../MafftAlignment/' + '_'.join(paralog) +'/' + '_'.join(paralog) +'_input.fasta'

            tree_newick = '../YeastTree.newick'
            DupLosList = '../YeastTestDupLost.txt'
            terminal_node_list = ['kluyveri', 'castellii', 'bayanus', 'kudriavzevii', 'mikatae', 'paradoxus', 'cerevisiae']
            node_to_pos = {'D1':0}
            pm_model = 'HKY'
            
            IGC_pm = 'One rate'

            if rate_variation:
                save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
                summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
                x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
            else:
                save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
                summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
                x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])


            test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                                 rate_variation, node_to_pos, terminal_node_list, save_file)
            test_JS.get_mle()
            #test_JS.get_expectedNumGeneconv()
            test_JS.get_individual_summary(summary_file)
            lnL_summary_file = summary_file.replace('summary.txt', 'lnL_summary.txt')
            test_JS.get_sitewise_loglikelihood_summary(lnL_summary_file)
