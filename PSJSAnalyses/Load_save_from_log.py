# This file is used to prepare necessary files for running PSJSGeneconv analysis
# Data: Yeast datasets
# Xiang Ji
# xji3@ncsu.edu
import os
import numpy as np

if __name__ == '__main__':
    pairs = []
    all_pairs = '../Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    for pair in pairs:
        for guess in [1, 2]:
            log_file = "./log/PSJS_HKY_" + "_".join(pair) + '_One_rate_Guess_' + str(guess) + '_rv_SCOK_nonclock_log.txt'
            save_file = './save/PSJS_HKY_' + "_".join(pair) + '_One_rate_Guess_' + str(guess) + '_rv_SCOK_nonclock_save.txt'

            with open(log_file, 'r') as f:
                all_log = f.readlines()[1:]
                num_items = set([len(line.split('\t')) for line in all_log])
                all_obj = [abs(float(line.split('\t')[0])) for line in all_log]
                best_line_num = all_obj.index(min(all_obj))

                with open(save_file, 'w+') as g:
                    g.write('\n'.join(all_log[best_line_num].split('\t')[1:min(num_items)]) + '\n')
            
