from .kmer import distance_matrix
from .tree import upgma
from .progressive import progressive

def align(sequences, names, k=3):
    matrix = distance_matrix(sequences, k)
    tree, merge_order = upgma(matrix, names)
    msa, msa_names = progressive(sequences, names, merge_order)
    return msa, msa_names