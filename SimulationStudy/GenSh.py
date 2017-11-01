import os

if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    sh_line = 'sbatch -o HMM-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Sim_'

    
    mean_tract_list = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]

##    sh_line = 'sbatch -o HKY-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/HKY_'
##    for geo in mean_tract_list:
##        IGC_bash_file = './HKY_Tract_' + str(geo) + '.sh'
##        with open(IGC_bash_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for sim_num in range(1, 101):
##                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
##                with open('./ShFiles/HKY_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python Run_HKY_PSJS.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity \n')

##    sh_line = 'sbatch -o HKY-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Old_'
##    for geo in mean_tract_list:
##        IGC_bash_file = './Old_Tract_' + str(geo) + '.sh'
##        with open(IGC_bash_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for sim_num in range(1, 101):
##                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
##                with open('./ShFiles/Old_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python Run_Old.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity \n')

##    sh_line = 'sbatch -o Tenth-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Tenth_'
##    for geo in mean_tract_list:
##        IGC_bash_file = './Tenth_Tract_' + str(geo) + '.sh'
##        with open(IGC_bash_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for sim_num in range(1, 101):
##                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
##                with open('./ShFiles/Tenth_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python Run_HKY_PSJS.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity --Case Tenth \n')
##
##    sh_line = 'sbatch -o Half-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Half_'
##    for geo in mean_tract_list:
##        IGC_bash_file = './Half_Tract_' + str(geo) + '.sh'
##        with open(IGC_bash_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for sim_num in range(1, 101):
##                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
##                with open('./ShFiles/Half_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python Run_HKY_PSJS.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity --Case Half \n')

##    sh_line = 'sbatch -o Bias-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Bias_HKY_'
##    for geo in mean_tract_list:
##        IGC_bash_file = './Bias_HKY_Tract_' + str(geo) + '.sh'
##        with open(IGC_bash_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for sim_num in range(1, 101):
##                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
##                with open('./ShFiles/Bias_HKY_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python Run_PSJS_Bias.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity \n')


##    sh_line = 'sbatch -o Grid-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Grid_HKY_'
##    for geo in mean_tract_list:
##        IGC_bash_file = './Grid_HKY_Tract_' + str(geo) + '.sh'
##        with open(IGC_bash_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for sim_num in range(1, 101):
##                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
##                with open('./ShFiles/Grid_HKY_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python Run_PSJS_Grid.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity \n')


    sh_line = 'sbatch -o True-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/True_'
    for geo in mean_tract_list:
        IGC_bash_file = './True_Tract_' + str(geo) + '.sh'
        with open(IGC_bash_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for sim_num in range(1, 101):
                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
                with open('./ShFiles/True_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    g.write('python Run_PSJS_TrueValue.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity \n')
