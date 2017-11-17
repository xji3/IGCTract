# This file prepares input files to run my program on Harpak data
# Xiang Ji
# xji3@ncsu.edu
import numpy as np
import os

if __name__ == '__main__':
    file_list_file = './missing_0_species_list.txt'
    file_list = np.loadtxt(file_list_file, dtype = str)

    species_list = ['Human', 'Chimp', 'Goril', 'Orang', 'Macaq', 'Mouse']

    if not os.path.isdir('./prepared_input/'):
        os.mkdir('./prepared_input')
        
    for seq_file in file_list:
        row_num = 0
        with open('./intronAlignmentSeparated/' + seq_file, 'r') as f:
            with open('./prepared_input/' + seq_file.replace('.', '_').replace('_pos_seq_formatted', '.fasta'), 'w+') as g:
                with open('./prepared_input/' + \
                          seq_file.replace('.pos.seq.formatted', '.seq.index')\
                          .replace('.', '_') + '.txt' , 'w+') as h:
                    for line in f:
                        item = line.replace('\n', '').split('\t')
                        if row_num == 0:
                            for i in range(len(item[1])):
                                if item[1][i] == '-' or item[1][i].upper() == 'N':
                                    print seq_file
                                else:
                                    h.write('\t'.join([str(i+1), '0', '0\n']))
                                
                        name = species_list[row_num/2] + '__Paralog' + str(row_num%2 + 1)#item[0].split('.')[0]
                        g.write('>' + name + '\n')
                        g.write(item[1] + '\n')
                        row_num += 1
            
                

    with open('./GeneToOrlg.txt', 'w+') as f:
        for species in species_list:
            if not species == 'Mouse':
                for i in range(1, 3):
                    f.write(species + '__Paralog' + str(i) + '\t' + str(i) + '\n')
            else:
                f.write(species + '__Paralog1\t0\n')
    
