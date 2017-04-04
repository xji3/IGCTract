import os
import subprocess
from Bio import Seq, SeqIO, AlignIO

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

def processAlignment(input_file):
    align = AlignIO.read(input_file,'fasta')
    i=0
    while(i<align.get_alignment_length()):
        if not align[:,i].find('-') == -1:
            if i==0:
                align = align[:,1:]
            elif i==(align.get_alignment_length()-1):
                align = align[:,:-1]
            else:
                align = align[:,:i]+align[:,(i+1):]
        else:
            i=i+1
    assert(align.get_alignment_length()%3==0)
    return align

def GapRemovedFasta(align, output_fasta):
    with open(output_fasta, 'w+') as f:
        for rec in align:
            f.write('>' + str(rec.id) + '\n')
            f.write(str(rec.seq)[:-3] + '\n')

    
if __name__ == '__main__':
    path = '/Users/xji3/Genconv/MafftAlignment/'
    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    pairs.remove(['YLR028C', 'YMR120C'])
    for pair in pairs:
        mkdir_cmd = ['mkdir', '_'.join(pair)]
        MAFFT_cmd = ['/usr/local/bin/mafft', '--auto',
                     path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA.fa', '>',
                     path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT.fa']
        cp_cmd = ['cp',
                  '../PairsAlignemt/' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa',
                  './' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa']
        subprocess.call(mkdir_cmd)
        subprocess.call(cp_cmd)
        
        translateDNAtoAA('./' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa',
                         './' + '_'.join(pair) + '/' + '_'.join(pair) + '_AA.fa')
        #os.system(' '.join(MAFFT_cmd))
        format_fasta(path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT.fa',
                     path + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT_formated.fa')
        translateAAAlignmentDNAAlignment('./' + '_'.join(pair) + '/' + '_'.join(pair) + '_AA_MAFFT_formated.fa',
                                         './' + '_'.join(pair) + '/' + '_'.join(pair) + '.fa',
                                         './' + '_'.join(pair) + '/' + '_'.join(pair) + '_MAFFT.fa')

        GapRemovedFasta(processAlignment('./' + '_'.join(pair) + '/' + '_'.join(pair) + '_MAFFT.fa'),
                        './' + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta')
                                         
                                         
    
