from Bio import SeqIO
import os
import subprocess

def load_fasta(fil):
    ''' Read a fasta file and return ids, seqs'''
    IDs, seqs = [], []
    for record in SeqIO.parse(fil, "fasta"):
        IDs.append(record.id)
        seqs.append(str(record.seq))
    return IDs, seqs

def write_fasta(names, seqs, outfile):
    ''' Write a fasta file '''
    with open(outfile, 'w') as f:
        for name, seq in zip(names, seqs):
            f.write(f">{name}\n{seq}\n")

def clean_seqs(seqs):
    return [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]

def downsample_msa(input_file, seq_id=0.95):
    '''Downsample sequences in MSA using mmseqs easy-cluster'''
    print("PATH in Python:", os.environ.get('PATH'))

    sp_command = ['mmseqs', 'easy-cluster', input_file, 'out2',
                'tmp2', '-c', str(seq_id), '--cov-mode', '1',
                '--threads', '4']
    
    subprocess.run(sp_command, env = os.environ, check=True)
    filt_IDs, filt_seqs = load_fasta('out2_rep_seq.fasta')
    return filt_IDs, filt_seqs

# NEXT: need to write something to remove gaps in first sequence
# 
# 
def remove_gap_cols(seqs):
    def get_char_locations(gapped_string):
        locs=[]
        for i, char in enumerate(gapped_string):
            if char != '-':
                locs.append(i)
        return locs

    query_char_locations = get_char_locations(seqs[0])

    new_seqs=[]
    for seq in seqs:
        new_seqs.append(''.join([char for i, char in enumerate(seq) if i in query_char_locations]))

    return new_seqs
