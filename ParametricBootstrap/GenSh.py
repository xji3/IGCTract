import os

if __name__ == '__main__':
    pairs = []
    all_pairs = '../Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    sh_line = 'sbatch -o pb-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'
    num_partition = 4
    for partition in range(num_partition):
        partition_start = 100 / num_partition * partition + 1
        partition_stop = 100 / num_partition * (partition + 1) + 1
        IGC_bash_file = './PSJSAnalyses_Partition_' + str(partition + 1) + '.sh'
        with open(IGC_bash_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for pair in pairs:
                single_bash_file_name = '_'.join(pair) + '_PSJSAnalyses_Partition_' + str(partition + 1) + '.sh'
                f.write(sh_line + single_bash_file_name +' \n')
                with open('./ShFiles/' + single_bash_file_name, 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    for sim_num in range(partition_start, partition_stop):             
                        g.write('python Run_HKY_PSJS.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --sim_num ' + str(sim_num) + ' --heterogeneity --case PSJSAnalyses \n')

        IGC_bash_file = './JSAnalyses_Partition_' + str(partition + 1) + '.sh'
        with open(IGC_bash_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for pair in pairs:
                single_bash_file_name = '_'.join(pair) + '_JSAnalyses_Partition_' + str(partition + 1) + '.sh'
                f.write(sh_line + single_bash_file_name +' \n')
                with open('./ShFiles/' + single_bash_file_name, 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    for sim_num in range(partition_start, partition_stop):             
                        g.write('python Run_HKY_PSJS.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --sim_num ' + str(sim_num) + ' --heterogeneity --case JSAnalyses \n')
