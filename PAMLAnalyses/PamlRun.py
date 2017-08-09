import os
import subprocess
from Bio import Seq, SeqIO, AlignIO
from Bio.Phylo.PAML import codeml, baseml
import numpy as np


def initialize(paralog, out_path = './output/', alignment_path = '../MafftAlignment/'):
    if not os.path.isdir(out_path + '_'.join(paralog)):
        os.mkdir(out_path + '_'.join(paralog))

    input_alignment = out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta' 
    if not os.path.isfile(input_alignment):
        subprocess.check_output(['cp', input_alignment.replace(out_path, alignment_path), input_alignment])

    input_tree = out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_tree.newick'
    old_paml_tree_path = '/Users/xji3/Genconv_Copy/NewClusterPackRun/NewPairsAlignment/'
    if not os.path.isfile(input_tree):
        subprocess.check_output(['cp', input_tree.replace(out_path, old_paml_tree_path), input_tree])


def run_paml(paralog, ctl_file, codeml_dir = '/Users/xji3/Downloads/paml4.8/bin/codeml'):
    codeml_cmd = [codeml_dir, ctl_file]
    os.chdir('/Users/xji3/GitFolders/Genconv/PAMLCheck/output/' + '_'.join(paralog) + '/')
    print(codeml_cmd)
    subprocess.check_output(codeml_cmd)
    os.chdir('/Users/xji3/GitFolders/Genconv/PAMLCheck')

if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    
    for paralog in pairs:
        initialize(paralog)
        ctl_file = '/Users/xji3/GitFolders/YeastIGCTract/PAMLAnalyses/output/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_codeml.ctl'
        run_paml(paralog, ctl_file)

    summary_mat = []
    finished_list = []
    label = ['MG94_codeml_tree_length', 'MG94_codeml_lnL']
    footer = ' '.join(label)


    #pairs = pairs[0:2]
    for pair in pairs:
        codeml_result = codeml.read('/Users/xji3/GitFolders/Genconv/PAMLCheck/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml_result.txt')
        summary_mat.append([codeml_result['NSsites'][0]['tree length'],
                            codeml_result['NSsites'][0]['lnL']])
        finished_list.append(pair)

    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
    np.savetxt(open('/Users/xji3/GitFolders/YeastIGCTract/PAMLAnalyses/paml_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
