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
    os.remove('out2_all_seqs.fasta')
    os.remove('out2_cluster.tsv')
    return filt_IDs, filt_seqs

from typing import List, Tuple
import sys
import re
ALIGNED_AA = set("ACDEFGHIKLMNPQRSTVWY-")  # allowed aligned chars (upper + gap)

def parse_sequence(seq: str) -> List[Tuple[str, str]]:
    """
    Parse a single A3M sequence string into a list of (aligned_char, insertions_after).
    - aligned_char: uppercase letter or '-' (gap)
    - insertions_after: string (possibly empty) of lowercase letters that follow this aligned char

    This version tolerates initial leading lowercase insertions by assigning them to
    a synthetic leading aligned gap '-' (so columns line up).
    """
    result = []
    i = 0
    # if sequence starts with lowercase insertions, attach them to a leading '-'
    if i < len(seq) and seq[i].islower():
        # consume them and store as insertions after a leading gap
        ins = []
        while i < len(seq) and seq[i].islower():
            ins.append(seq[i])
            i += 1
        result.append(('-', ''.join(ins)))

    while i < len(seq):
        ch = seq[i]
        if ch.isupper() or ch == '-':
            aligned_char = ch
            i += 1
            # collect following lowercase insertions
            ins = []
            while i < len(seq) and seq[i].islower():
                ins.append(seq[i])
                i += 1
            result.append((aligned_char, ''.join(ins)))
        else:
            # unexpected case (shouldn't happen now), but handle gracefully
            raise ValueError(f"Unexpected lowercase at position {i} in sequence: {seq!r}")
    return result


def expand_sequences(parsed_seqs: List[List[Tuple[str, str]]]) -> List[str]:
    """
    Given parsed sequences (list of parse_sequence outputs, one per sequence),
    compute per-aligned-position maximum insertion length, and expand each sequence
    into the full column-wise representation.

    Returns list of expanded sequences (strings). Inserted letters are kept (lowercase).
    For sequences lacking insertion at a given slot, '-' are inserted to pad to the max length.
    """
    # sanity: all sequences must have the same number of aligned positions
    lengths = [len(p) for p in parsed_seqs]
    if len(set(lengths)) != 1:
        raise ValueError(f"Parsed sequences have differing numbers of aligned positions: {lengths}")

    npos = lengths[0]
    # compute max insertion length after each aligned position
    max_ins = [0] * npos
    for p in parsed_seqs:
        for i, (_, ins) in enumerate(p):
            if len(ins) > max_ins[i]:
                max_ins[i] = len(ins)

    # build expanded sequences
    expanded = []
    for p in parsed_seqs:
        out_chars = []
        for i, (aligned_char, ins) in enumerate(p):
            # put the aligned char (keep '-' or uppercase as-is)
            out_chars.append(aligned_char)
            # append insertion letters (if any), then pad with '-' to max_ins[i]
            if ins:
                out_chars.extend(list(ins.upper()))
                pad = max_ins[i] - len(ins)
                if pad > 0:
                    out_chars.extend(['-'] * pad)
            else:
                # no insertion in this seq at this position -> add all gaps
                if max_ins[i] > 0:
                    out_chars.extend(['-'] * max_ins[i])
        expanded.append(''.join(out_chars))
    return expanded

def remove_first_gaps(seqs):
    """
    Remove columns where seqs[0] has '-' and convert other sequences' residues
    in those columns into lowercase insertions attached to the previous retained column.

    Returns (headers, new_seqs).
    """
    if len(seqs) == 0:
        return headers, seqs

    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise ValueError("All sequences must be the same length")

    n = len(seqs)
    first = seqs[0]

    outs: List[List[str]] = [[] for _ in range(n)]
    ins_buffers: List[List[str]] = [[] for _ in range(n)]

    any_retained = False  # have we encountered any retained aligned column yet?

    for col in range(L):
        if first[col] != '-':
            # retained aligned column: first has an aligned char
            any_retained = True
            for j in range(n):
                if ins_buffers[j]:
                    # append the collected insertions (already lowercased)
                    outs[j].extend(ins_buffers[j])
                    ins_buffers[j] = []
                # append this aligned character (keep as-is, so uppercase or '-')
                outs[j].append(seqs[j][col])
        else:
            for j in range(1, n):  # skip j==0 (first sequence) since it's the one driving removal
                ch = seqs[j][col]
                if ch != '-':
                    # convert to lowercase and stash into buffer to be appended after previous retained
                    ins_buffers[j].append(ch.lower())

    for j in range(n):
        if ins_buffers[j]:
            outs[j].extend(ins_buffers[j])
            ins_buffers[j] = []

    new_seqs = [''.join(lst) for lst in outs]

    return new_seqs

def fix_neighborcluster_msas(seqs):
    '''seqs: list of strings'''
    parsed = [parse_sequence(s) for s in seqs]
    expanded = expand_sequences(parsed)
    final = remove_first_gaps(expanded)
    return final