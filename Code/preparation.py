# This file is used to prepare necessary files for running PSJSGeneconv analysis
# Data: Yeast datasets
# Xiang Ji
# xji3@ncsu.edu
import os

if __name__ == '__main__':
 
    pairs = []
    all_pairs = '../All_Pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    pairs.remove(['YLR028C', 'YMR120C'])

    # prepare GeneToOrlg files
    if not os.path.isdir('../GeneToOrlg'):
        os.mkdir('../GeneToOrlg')

    ingroup_node_list = ['cerevisiae', 'paradoxus', 'mikatae', 'kudriavzevii', 'bayanus', 'castellii']
    for paralog in pairs:
        with open('../GeneToOrlg/' + '_'.join(paralog) + '_GeneToOrlg.txt', 'w+') as f:
            for node in ingroup_node_list:
                f.write('__'.join([node, paralog[0]]) + '\t 1 \n')
                f.write('__'.join([node, paralog[1]]) + '\t 2 \n')

            f.write('__'.join(['kluyveri', paralog[0]]) + '\t 0 \n')

    # Gen sh files
    if not os.path.isdir('../ShFiles'):
        os.mkdir('../ShFiles')

    IGC_pm = 'One_rate'
    sh_line = 'sbatch -o PSJS-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'
    sh_file_all_name = './PSJS_' + IGC_pm + '_all.sh'
    with open(sh_file_all_name, 'w+') as g:
        g.write('#!/bin/bash' + '\n')
        for paralog in pairs:
            sh_file_name = 'PSJS_HKY_' + '_'.join(paralog) + '_' + IGC_pm + '_nonclock.sh'
            with open('../ShFiles/' + sh_file_name, 'w+') as f:
                f.write('#!/bin/bash' + '\n')
                f.write('python Run.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + '\n')
            g.write(sh_line + sh_file_name + '  \n')
                
            
