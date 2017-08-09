import os
import subprocess
from Bio import Seq, SeqIO, AlignIO
from Bio.Phylo.PAML import codeml, baseml
import numpy as np

def partition_seqfile(seqfile, partitioned_seqfile):
    assert(os.path.isfile(seqfile))
    align = AlignIO.read(seqfile, 'fasta')
    align_length = align.get_alignment_length()
    with open(partitioned_seqfile, 'w+') as f:
        f.write('    '.join([' ', '13', str(align_length), 'GC \n']))
        with open(seqfile, 'r') as g:
            for line in g:
                if line[0] == '>':
                    name = line[1:-1]
                else:
                    f.write(name + '    ' + line)
                
    

if __name__ == '__main__':
    path = '/Users/xji3/GitFolders/YeastIGCTract/PAMLAnalyses/'
    
    pair_file = './Filtered_pairs.txt'
    pairs = []
    with open(pair_file, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    tree_pair = ['YML026C', 'YDR450W']
    with open(path + 'YeastTree.newick', 'r') as f:
        all_tree_lines = f.readlines()

    with open(path + 'codeml_tail.ctl', 'r') as f:
        all_codeml_ctl_lines = f.readlines()

    with open(path + 'baseml_tail.ctl', 'r') as f:
        all_baseml_ctl_lines = f.readlines()
    

    codeml_dir = '/Users/xji3/Downloads/paml4.8/bin/codeml'
    baseml_dir = '/Users/xji3/Downloads/paml4.8/bin/baseml'

    #pairs = [pairs[0]]
    for pair in pairs:
        print 'Now run paml on pair ' + ' '.join(pair)
        seqfile = path + '/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta'
        partitioned_seqfile = seqfile.replace('input.fasta', 'partitioned.fasta')
        if not os.path.isdir(path + 'output/' + '_'.join(pair)):
            os.mkdir(path + 'output/' + '_'.join(pair))
            
        if not os.path.isfile(seqfile):
            mafft_aligned_file = '../MafftAlignment/' + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta'
            os.system(' '.join(['cp', mafft_aligned_file, seqfile]))

        #if not os.path.isfile(partitioned_seqfile):
        partition_seqfile(seqfile, partitioned_seqfile)
            
        treefile = path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_tree.newick'        
        with open(treefile, 'w+') as f:
            for line in all_tree_lines:
                new_line = line.replace(tree_pair[0], '__'+pair[0])
                new_line = new_line.replace(tree_pair[1], '__'+pair[1])
                f.write(new_line)


        outfile_baseml = path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml'
        baseml_ctlfile = path + 'output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml_control.ctl'
        with open(baseml_ctlfile, 'w+') as f:
            f.writelines(['seqfile = ' + partitioned_seqfile + '\n', 'treefile = ' + treefile + '\n', 'outfile = ' + outfile_baseml + '\n'])
            f.writelines(all_baseml_ctl_lines)

##        baseml_cmd = [baseml_dir, './output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml_control.ctl']
##        subprocess.check_output(baseml_cmd)
        

    summary_mat = []
    finished_list = []
    label = ['HKY_baseml_tree_length', 'HKY_baseml_lnL', 'HKY_baseml_kappa', 'HKY_r2', 'HKY_r3']
    footer = ' '.join(label)


    #pairs = pairs[0:2]
    for pair in pairs:
        #codeml_result = codeml.read('/Users/xji3/Genconv_Copy/NewClusterPackRun/NewPairsAlignment/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml')
        baseml_result = baseml.read('/Users/xji3/GitFolders/YeastIGCTract/PAMLAnalyses/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
        parameter_list = baseml_result['parameters']['parameter list'].split(' ')
        summary_mat.append([baseml_result['tree length'],
                            baseml_result['lnL'],
                            float(parameter_list[-1]),
                            float(parameter_list[-3]),
                            float(parameter_list[-2])]
                           )
        finished_list.append(pair)

    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
    np.savetxt(open('/Users/xji3/GitFolders/YeastIGCTract/PAMLAnalyses/output/paml_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
