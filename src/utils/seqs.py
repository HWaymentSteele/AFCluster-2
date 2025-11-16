import numpy as np

def encode_seqs(seqs, max_len=108, alphabet=None):
    if alphabet is None:
        alphabet = "ACDEFGHIKLMNPQRSTVWY-"
    A = len(alphabet)
    lookup = {res: i for i, res in enumerate(alphabet)}

    arr = np.zeros((len(seqs), max_len, A), dtype=np.float32)
    for j, seq in enumerate(seqs):
        idxs = [lookup.get(c, A-1) for c in seq[:max_len]]  # unknowns -> '-'
        arr[j, np.arange(len(idxs)), idxs] = 1
    return arr.reshape(len(seqs), max_len * A)

def get_closest_n_seqs(seq, other_seqs, n=30):
    ohe_seq = encode_seqs([seq], max_len=len(seq))
    ohe_other_seqs = encode_seqs(other_seqs, max_len=len(seq))

    dists = np.linalg.norm(ohe_other_seqs - ohe_seq, axis=1)
    closest_indices = np.argsort(dists)[:n+1] #note this is also including the original seq, that's fine
    
    return closest_indices