import os

if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    sh_line = 'sbatch -o HMM-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Sim_'

    
    mean_tract_list = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]
    for geo in mean_tract_list:
        IGC_bash_file = './Sim_Tract_' + str(geo) + '.sh'
        with open(IGC_bash_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for sim_num in range(1, 101):
                f.write(sh_line + str(sim_num) + '_Tract_' + str(geo) + '.sh \n')
                with open('./ShFiles/Sim_' + str(sim_num) + '_Tract_' + str(geo) + '.sh', 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    g.write('python Run.py --geo ' + str(geo) + ' --sim_num ' + str(sim_num) + ' --heterogeneity \n')


##    sh_line = 'sbatch -o IGCSim-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles_multiply3/'
##
##    for IGC_geo in IGC_geo_list:
##        IGC_geo_sh_file = './' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_multiply3.sh'
##        with open(IGC_geo_sh_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for num_sim in sim_num_list:
##                f.write(sh_line + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh \n')
##                with open('./ShFiles_multiply3/' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python RunSimulation_multiply3.py --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')
