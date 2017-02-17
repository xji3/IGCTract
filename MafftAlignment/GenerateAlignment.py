import os
import subprocess
from Bio import Seq, SeqIO, AlignIO
import numpy as np
from math import floor

def translateDNAtoAA(input_fasta, output_fasta):  
    with open(input_fasta, 'r') as f:
        with open(output_fasta, 'w+') as g:
            for line in f.readlines():
                if line[0] == '>':
                    g.write(line)
                    continue
                else:
                    assert(len(line) %3 == 1)
                    g.write(Seq.translate(line[:-1], to_stop = True) + '\n')

def format_fasta(input_fasta, output_fasta):
    with open(input_fasta, 'r') as f:
        with open(output_fasta, 'w+') as g:
            seq = ''
            for line in f.readlines():
                if line[0] == '>' and seq == '': # first line
                    g.write(line)
                elif line[0] == '>' and seq != '':
                    g.write(seq + '\n')
                    g.write(line)
                    seq = ''
                else:
                    seq += line[:-1]
            g.write(seq + '\n')
                    
    
def translateAAAlignmentDNAAlignment(AA_alignment, DNA_fasta, output_fasta):
    seq_dict = SeqIO.to_dict(SeqIO.parse( DNA_fasta, "fasta" ))
    name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
    with open(AA_alignment, 'r') as f:
        with open(output_fasta, 'w+') as g:
            for line in f.readlines():
                if line[0] == '>':
                    name = line[1:-1]
                    g.write(line)
                else:
                    dna_seq = name_to_seq[name]
                    gap = 0
                    new_line = ''
                    for i in range(len(line) -1):
                        if line[i] == '-':
                            new_line += '---'
                            gap += 1
                        else:
                            new_line += dna_seq[(3 * (i - gap)):(3 * (i - gap) + 3)]
                    g.write(new_line + '\n')

def processAlignment(input_file, reference_seq_name):
    align = AlignIO.read(input_file,'fasta')
    # now get reference_seq_row_num
    ref_row_num = [rec.id for rec in align].index(reference_seq_name)    
    i=0
    ref_pos = 0
    seq_index = []
    while(i<align.get_alignment_length()):
        if not align[:,i].find('-') == -1:
            if not align[ref_row_num, i] == '-':
                # ref seq does not have gap in this column
                # this column is deleted and ref_pos should move by 1
                ref_pos += 1
                
            if i==0:
                align = align[:,1:]
            elif i==(align.get_alignment_length()-1):
                align = align[:,:-1]
            else:
                align = align[:,:i]+align[:,(i+1):]
        else:
            i=i+1
            ref_pos += 1
            seq_index.append([ref_pos, floor((ref_pos - 0.5) / 3.0) + 1, ref_pos - floor((ref_pos - 0.5) / 3.0) * 3.0])
    assert(align.get_alignment_length()%3==0)
    return align, seq_index

def GapRemovedFasta(align, output_fasta):
    with open(output_fasta, 'w+') as f:
        for rec in align:
            spe = str(rec.id)[:-7]
            paralog = str(rec.id)[-7:]
            f.write('>' + spe + '__' + paralog + '\n')
            f.write(str(rec.seq)[:-3] + '\n')

def get_seq_index(name_to_seq, ref_seq_name, ref_input_file, idx_seq_file):
    ref_seq_dict = SeqIO.to_dict(SeqIO.parse(ref_input_file, 'fasta'))
    ref_name_to_seq = {name:str(ref_seq_dict[name].seq) for name in ref_seq_dict.keys()}
    ref_seq = ref_name_to_seq[ref_seq_name]
    
    idx_seq_aligned_file = idx_seq_file.replace('ref_seq', 'ref_seq_aligned')
    with open(idx_seq_file, 'w+') as f:
        f.write('>' + ref_seq_name + '\n')
        f.write(ref_seq + '\n')
        f.write('>' + ref_seq_name + '_ref\n')
        f.write(name_to_seq[pair[0]] + '\n')

    # now perform alignment using MAFFT
    mafft_cmd = ['/usr/local/bin/mafft', '--auto', idx_seq_file, '>', idx_seq_aligned_file]
    #os.system(' '.join(mafft_cmd))

    # now generate seq_idx file
    align = AlignIO.read(idx_seq_aligned_file,'fasta')
    # now get reference_seq_row_num
    ref_row_num = [rec.id for rec in align].index(ref_seq_name)    
    i=0
    ref_pos = 0
    seq_index = []
    while(i<align.get_alignment_length()):
        if not align[:,i].find('-') == -1:
            if align[ref_row_num, i] == '-':
                # ref seq does not have gap in this column
                # this column is deleted and ref_pos should move by 1
                ref_pos += 1                
            if i==0:
                align = align[:,1:]
            elif i==(align.get_alignment_length()-1):
                align = align[:,:-1]
            else:
                align = align[:,:i]+align[:,(i+1):]
        else:
            i=i+1
            ref_pos += 1
            seq_index.append([ref_pos, int(floor((ref_pos - 0.5) / 3.0) + 1), int(ref_pos - floor((ref_pos - 0.5) / 3.0) * 3.0)])
    return seq_index

    
if __name__ == '__main__':
    path = '/Users/xji3/GitFolders/YeastIGCTract/MafftAlignment/'
    #path = '/Users/Xiang/GitFolders/YeastIGCTract/MafftAlignment/'
    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    pairs.remove(['YLR028C', 'YMR120C'])

    dna_seq_file = '/Users/xji3/Documents/YeastGenome/orf_genomic.fasta'
    dna_seq_dict = SeqIO.to_dict(SeqIO.parse(dna_seq_file, 'fasta'))
    dna_name_to_seq = {name:str(dna_seq_dict[name].seq) for name in dna_seq_dict.keys()}

    cdna_seq_file = '/Users/xji3/Documents/YeastGenome/orf_coding.fasta'
    cdna_seq_dict = SeqIO.to_dict(SeqIO.parse(cdna_seq_file, 'fasta'))
    cdna_name_to_seq = {name:str(cdna_seq_dict[name].seq) for name in cdna_seq_dict.keys()}

    pair = ['YLR406C', 'YDL075W']
    ref_input_file = './' + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta'
    ref_seq_name = 'cerevisiae__' + pair[0]

    idx_seq_file = path + '_'.join(pair) + '/' + pair[0] + '_ref_seq.fasta'

    #seq_index = get_seq_index(name_to_seq, ref_seq_name, ref_input_file, idx_seq_file)

    ref_seq_dict = SeqIO.to_dict(SeqIO.parse(ref_input_file, 'fasta'))
    ref_name_to_seq = {name:str(ref_seq_dict[name].seq) for name in ref_seq_dict.keys()}
    ref_seq = ref_name_to_seq[ref_seq_name]
    
    idx_seq_aligned_file = idx_seq_file.replace('ref_seq', 'ref_seq_aligned')
    with open(idx_seq_file, 'w+') as f:
#        f.write('>' + ref_seq_name + '\n')
#        f.write(ref_seq + '\n')
        f.write('>' + ref_seq_name + '_dna\n')
        f.write(dna_name_to_seq[pair[0]] + '\n')
        f.write('>' + ref_seq_name + '_cdna\n')
        f.write(cdna_name_to_seq[pair[0]] + '\n')

    # now perform alignment using MAFFT
    mafft_cmd = ['/usr/local/bin/mafft', '--auto', idx_seq_file, '>', idx_seq_aligned_file]
    os.system(' '.join(mafft_cmd))

    # now generate seq_idx file
    align = AlignIO.read(idx_seq_aligned_file,'fasta')
    # now get reference_seq_row_num
    ref_row_num = [rec.id for rec in align].index(ref_seq_name)    
    i=0
    ref_pos = 0
    seq_index = []
    while(i<align.get_alignment_length()):
        if not align[:,i].find('-') == -1:
            if align[ref_row_num, i] == '-':
                # ref seq does not have gap in this column
                # this column is deleted and ref_pos should move by 1
                ref_pos += 1                
            if i==0:
                align = align[:,1:]
            elif i==(align.get_alignment_length()-1):
                align = align[:,:-1]
            else:
                align = align[:,:i]+align[:,(i+1):]
        else:
            i=i+1
            ref_pos += 1
            seq_index.append([ref_pos, int(floor((ref_pos - 0.5) / 3.0) + 1), int(ref_pos - floor((ref_pos - 0.5) / 3.0) * 3.0)])
   

##    for pair in pairs:
##        mkdir_cmd = ['mkdir', '_'.join(pair)]
##        MAFFT_cmd = ['/usr/local/bin/mafft', '--auto',
##                     path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA.fa', '>',
##                     path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT.fa']
##        cp_cmd = ['cp',
##                  '../PairsAlignemt/' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa',
##                  './' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa']
##        subprocess.call(mkdir_cmd)
##        subprocess.call(cp_cmd)
##        
##        translateDNAtoAA('./' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa',
##                         './' + '_'.join(pair) + '/' + '_'.join(pair) + '_AA.fa')
##        #os.system(' '.join(MAFFT_cmd))
##        format_fasta(path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT.fa',
##                     path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT_formated.fa')
##        translateAAAlignmentDNAAlignment('./' + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT_formated.fa',
##                                         './' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa',
##                                         './' + '_'.join(pair) + '/' + '_'.join(pair) + '_MAFFT.fa')
##
##        reference_seq_name = 'cerevisiae' + pair[0]
##        processed_align, seq_index = processAlignment('./' + '_'.join(pair) + '/' + '_'.join(pair) + '_MAFFT.fa', reference_seq_name)
##        np.savetxt('./' + '_'.join(pair) + '/' + '_'.join(pair) + '_seq_index.txt', np.array(seq_index[:-3], dtype = int), fmt = '%i')
##        GapRemovedFasta(processed_align, './' + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta')
                                         
                                         
##    reference_seq_name = 'cerevisiae' + pair[0]
##    reference_seq_name = 'kluyveri' + pair[0]
##    input_file = './' + '_'.join(pair) + '/' + '_'.join(pair) + '_MAFFT.fa'
##    align = AlignIO.read(input_file,'fasta')
##    # now get reference_seq_row_num
##    ref_row_num = [rec.id for rec in align].index(reference_seq_name)
##    
##    i=0
##    ref_pos = 0
##    seq_index = []
##    while(i<align.get_alignment_length()):
##        if not align[:,i].find('-') == -1:
##            if not align[ref_row_num, i] == '-':
##                # ref seq does not have gap in this column
##                # this column is deleted and ref_pos should move by 1
##                ref_pos += 1
##                
##            if i==0:
##                align = align[:,1:]
##            elif i==(align.get_alignment_length()-1):
##                align = align[:,:-1]
##            else:
##                align = align[:,:i]+align[:,(i+1):]
##        else:
##            i=i+1
##            ref_pos += 1
##            seq_index.append(ref_pos)
##    assert(align.get_alignment_length()%3==0)
