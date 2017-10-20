import os

if __name__ == '__main__':
    for geo in [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]:
        for sim_num in range(1,101):
            original_alignment = './Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) \
                                 + '/YDR418W_YEL054C_sim_' + str(sim_num) + '_newformat.fasta'
            half_alignment_first = original_alignment.replace('.fasta', '_first_half.fasta')
            half_alignment_second = original_alignment.replace('.fasta', '_second_half.fasta')
            with open(original_alignment, 'r') as f:
                with open(half_alignment_first, 'w+') as g1:
                    with open(half_alignment_second, 'w+') as g2:
                        for line in f:
                            if line[0] == '>':
                                g1.write(line)
                                g2.write(line)
                            else:
                                content = line.replace('\n', '')
                                assert(len(content) == 489)
                                g1.write(content[:243] + '\n')
                                g2.write(content[243:] + '\n')
                                
                    
