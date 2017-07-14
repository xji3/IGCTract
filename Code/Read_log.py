# This file is to change current parameter location according to log file's lowest point
# Xiang Ji
# xji3@ncsu.edu
import os
import numpy as np

def read_log(log_file):
    assert(os.path.isfile(log_file))
    logs = np.loadtxt(log_file)
    best_line = logs[:, 0].argmin()
    total_len = logs.shape[1]
    parameter_values = logs[best_line, 1:(total_len + 1)/2]
    return parameter_values

def overwrite_save_file(save_file, parameter_values):
    np.savetxt(open(save_file, 'w+'), parameter_values.T)
    
    

if __name__ == '__main__':
    pairs = []
    all_pairs = '../Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

 
    for pair in pairs[:]:
        for guess_num in range(1,3):
            log_file = './log/PSJS_HKY_' + '_'.join(pair) + '_One_rate_Guess_' + str(guess_num) + '_nonclock_log.txt'
            save_file = './save/PSJS_HKY_' + '_'.join(pair) + '_One_rate_Guess_' + str(guess_num) + '_nonclock_save.txt'
            parameter_values = read_log(log_file)
            overwrite_save_file(save_file, parameter_values)

            log_file = './log/PSJS_HKY_' + '_'.join(pair) + '_One_rate_Guess_' + str(guess_num) + '_rv_SCOK_nonclock_log.txt'
            save_file = './save/PSJS_HKY_' + '_'.join(pair) + '_One_rate_Guess_' + str(guess_num) + '_rv_SCOK_nonclock_save.txt'
            parameter_values = read_log(log_file)
            overwrite_save_file(save_file, parameter_values)

            # copy save file over
            cp_cmd = ['cp', save_file, save_file.replace('SCOK', 'NOSC')]
            
