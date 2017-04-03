# This file is used to prepare necessary files for running PSJSGeneconv analysis
# Data: Yeast datasets
# Xiang Ji
# xji3@ncsu.edu
import os
import numpy as np

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

    pairs = []
    all_pairs = '../Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    # Now get an averaged parameter value list for initial guess
    averaged_parameters = 0.0
    for pair in pairs:
        JS_save_file = './save/JS_HKY_' + '_'.join(pair) + '_One_rate_nonclock_save.txt'
        single_parameters = np.exp(np.loadtxt(JS_save_file))
        averaged_parameters += single_parameters
    averaged_parameters = averaged_parameters / len(pairs)
    np.savetxt('./averaged_JS_HKY_One_rate_nonclock_save.txt', np.log(averaged_parameters))
    

    IGC_pm = 'One_rate'
    
    sh_line = 'sbatch -p long -o PSJS-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ../ShFiles/'
    for tract_length in [ 5.0, 30.0, 200.0]:
        sh_file_all_name = './PSJS_' + IGC_pm +'_init_' + str(tract_length) + '_all.sh'
        with open(sh_file_all_name, 'w+') as g:
            g.write('#!/bin/bash' + '\n')
            for paralog in pairs:
                sh_file_name = '_'.join(paralog) + '_PSJS_HKY_' + IGC_pm +'_init_' + str(tract_length) +  '_nonclock.sh'
                with open('../ShFiles/' + sh_file_name, 'w+') as f:
                    f.write('#!/bin/bash' + '\n')
                    f.write('python Run.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --L ' + str(tract_length) + '\n')
                g.write(sh_line + sh_file_name + '  \n')

    sh_line = 'sbatch -p long -o PSJSG-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ../ShFiles/'
    for guess in [ 1, 2]:
        sh_file_all_name = './PSJS_' + IGC_pm +'_guess_' + str(guess) + '_all.sh'
        with open(sh_file_all_name, 'w+') as g:
            g.write('#!/bin/bash' + '\n')
            for paralog in pairs:
                sh_file_name = '_'.join(paralog) + '_PSJS_HKY_' + IGC_pm +'_guess_' + str(guess) +  '_nonclock.sh'
                with open('../ShFiles/' + sh_file_name, 'w+') as f:
                    f.write('#!/bin/bash' + '\n')
                    f.write('python Guess.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --G ' + str(guess)\
                            + ' --homogeneity --coding --samecodon \n')
                g.write(sh_line + sh_file_name + '  \n')

    sh_line = 'sbatch -p long -o PSJSG-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ../ShFiles/'
    for guess in [ 1, 2]:
        sh_file_all_name = './PSJS_' + IGC_pm +'_RV_guess_' + str(guess) + '_all.sh'
        with open(sh_file_all_name, 'w+') as g:
            g.write('#!/bin/bash' + '\n')
            for paralog in pairs:
                sh_file_name = '_'.join(paralog) + '_PSJS_HKY_' + IGC_pm +'_RV_guess_' + str(guess) +  '_nonclock.sh'
                with open('../ShFiles/' + sh_file_name, 'w+') as f:
                    f.write('#!/bin/bash' + '\n')
                    f.write('python Guess.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --G ' + str(guess)\
                            + ' --heterogeneity --coding --no-samecodon \n')
                    f.write('python Guess.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --G ' + str(guess)\
                            + ' --heterogeneity --coding --samecodon \n')
                g.write(sh_line + sh_file_name + '  \n')


    IGC_pm = 'One_rate'
    tract_length = 30.0
    sh_line = 'sbatch -p long -o PSJS-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ../ShFiles/'
    for dim in [ 1, 2]:
        sh_file_all_name = './PSJS_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) + '_all.sh'
        with open(sh_file_all_name, 'w+') as g:
            g.write('#!/bin/bash' + '\n')
            for paralog in pairs:
                sh_file_name = '_'.join(paralog) + '_PSJS_HKY_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) +  '_nonclock.sh'
                with open('../ShFiles/' + sh_file_name, 'w+') as f:
                    f.write('#!/bin/bash' + '\n')
                    f.write('python Run.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --D ' + str(dim)\
                            + ' --homogeneity --coding --samecodon \n')
                g.write(sh_line + sh_file_name + '  \n')

        plot_sh_file_all_name = './Plot_PSJS_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) + '_all.sh'
        with open(plot_sh_file_all_name, 'w+') as g:
            g.write('#!/bin/bash' + '\n')
            for paralog in pairs:
                sh_file_name = 'plot_' + '_'.join(paralog) + '_PSJS_HKY_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) +  '_nonclock.sh'
                with open('../ShFiles/' + sh_file_name, 'w+') as f:
                    f.write('#!/bin/bash' + '\n')
                    f.write('python plot.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --homogeneity --coding --zoom --D ' + str(dim) + '\n')
                    f.write('python plot.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --homogeneity --coding --no-zoom --D ' + str(dim) + '\n')
                g.write(sh_line + sh_file_name + '  \n')            

    # Rate heterogeneity bash file
    IGC_pm = 'One_rate'
    tract_length = 30.0
    sh_line = 'sbatch -p long -o PSJS-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ../ShFiles/'
    for allow_same_codon in [True, False]:
        for dim in [ 1, 2]:
            if allow_same_codon:
                sh_file_all_name = './PSJS_RV_SCOK_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) + '_all.sh'
                plot_sh_file_all_name = './Plot_PSJS_RV_SCOK_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) + '_all.sh'
            else:
                sh_file_all_name = './PSJS_RV_NOSC_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) + '_all.sh'
                plot_sh_file_all_name = './Plot_PSJS_RV_NOSC_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) + '_all.sh'
                
            with open(sh_file_all_name, 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                for paralog in pairs:
                    if allow_same_codon:
                        sh_file_name = '_'.join(paralog) + '_PSJS_RV_SCOK_HKY_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) +  '_nonclock.sh'
                    else:
                        sh_file_name = '_'.join(paralog) + '_PSJS_RV_NOSC_HKY_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) +  '_nonclock.sh'
                        
                    with open('../ShFiles/' + sh_file_name, 'w+') as f:
                        f.write('#!/bin/bash' + '\n')
                        if allow_same_codon:
                            f.write('python Run.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --heterogeneity --coding --samecodon --D ' + str(dim) + '\n')
                        else:
                            f.write('python Run.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --heterogeneity --coding --no-samecodon --D ' + str(dim) + '\n')
                    g.write(sh_line + sh_file_name + '  \n')

            with open(plot_sh_file_all_name, 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                for paralog in pairs:
                    if allow_same_codon:
                        sh_file_name = 'plot_' + '_'.join(paralog) + '_PSJS_RV_SCOK_HKY_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) +  '_nonclock.sh'
                    else:
                        sh_file_name = 'plot_' + '_'.join(paralog) + '_PSJS_RV_NOSC_HKY_dim_' + str(dim) + '_' + IGC_pm +'_init_' + str(tract_length) +  '_nonclock.sh'
                        
                    with open('../ShFiles/' + sh_file_name, 'w+') as f:
                        f.write('#!/bin/bash' + '\n')
                        if allow_same_codon:
                            f.write('python plot.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --heterogeneity --coding --samecodon --zoom --D ' + str(dim) + '\n')
                            f.write('python plot.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --heterogeneity --coding --samecodon --no-zoom --D ' + str(dim) + '\n')
                        else:
                            f.write('python plot.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --heterogeneity --coding --no-samecodon --zoom --D ' + str(dim) + '\n')
                            f.write('python plot.py --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --heterogeneity --coding --no-samecodon --no-zoom --D ' + str(dim) + '\n')
                    g.write(sh_line + sh_file_name + '  \n')
