from .hmm import build_profile, viterbi

def insert_gaps(seq, gaps_in_prof):
    result = []
    seq_i = 0
    for gap in gaps_in_prof:
        if gap:
            result.append('-')  
        else:
            if seq_i < len(seq):
                result.append(seq[seq_i])
                seq_i += 1
    return ''.join(result)

def align_two(seq1, seq2):
    """İki diziyi hizalar"""
    profile = build_profile([seq1])
    seq2_aligned, gaps_in_prof = viterbi(profile, seq2)
    seq1_aligned = insert_gaps(seq1, gaps_in_prof)
    return seq1_aligned, seq2_aligned

def align_to_profile(aligned_seqs, new_seq):
    """Profil + yeni dizi hizala, mevcut dizileri güncelle"""
    profile = build_profile(aligned_seqs)
    new_aligned, gaps_in_prof = viterbi(profile, new_seq)

    updated = []
    for seq in aligned_seqs:
        updated.append(insert_gaps(seq, gaps_in_prof))

    return updated, new_aligned

def progressive(sequences, names, merge_order):
    seq_dict = {name: seq for name, seq in zip(names, sequences)}
    aligned_seqs_dict  = {}
    aligned_names_dict = {}

    print("\n=== Progressive Alignment ===")

    for pair in merge_order:
        label1, label2 = pair

        if label1 in aligned_seqs_dict:
            seqs1 = aligned_seqs_dict[label1]
            nams1 = aligned_names_dict[label1]
        elif label1 in seq_dict:
            seqs1 = [seq_dict[label1]]
            nams1 = [label1]
        else:
            seqs1 = []
            nams1 = []

        if label2 in aligned_seqs_dict:
            seqs2 = aligned_seqs_dict[label2]
            nams2 = aligned_names_dict[label2]
        elif label2 in seq_dict:
            seqs2 = [seq_dict[label2]]
            nams2 = [label2]
        else:
            seqs2 = []
            nams2 = []

        print(f"Hizalanıyor: {label1} + {label2}")

        if len(seqs1) == 1 and len(seqs2) == 1:
            s1, s2 = align_two(seqs1[0], seqs2[0])
            result_seqs  = [s1, s2]
            result_names = [nams1[0], nams2[0]]
        elif len(seqs1) == 1 and len(seqs2) > 1:
            updated, new_aligned = align_to_profile(seqs2, seqs1[0])
            result_seqs  = [new_aligned] + updated
            result_names = [nams1[0]] + list(nams2)
        else:
            result_seqs  = list(seqs1)
            result_names = list(nams1)
            for seq, name in zip(seqs2, nams2):
                updated, new_aligned = align_to_profile(result_seqs, seq)
                result_seqs  = updated + [new_aligned]
                result_names = result_names + [name]

        new_label = f"({label1},{label2})"
        aligned_seqs_dict[new_label]  = result_seqs
        aligned_names_dict[new_label] = result_names

    final_label = list(aligned_seqs_dict.keys())[-1]
    final_seqs  = aligned_seqs_dict[final_label]
    final_names = aligned_names_dict[final_label]

    print("\n=== Final MSA ===")
    for name, seq in zip(final_names, final_seqs):
        print(f"{name}: {seq}")

    return final_seqs, final_names


if __name__ == "__main__":
    from kmer import distance_matrix
    from tree import upgma

    sequences = ["ACGTAGCT", "ACTAGCT", "TTGCAGCT"]
    names     = ["Dizi1", "Dizi2", "Dizi3"]

    matrix = distance_matrix(sequences)
    tree, merge_order = upgma(matrix, names)
    print(f"Merge sirasi: {merge_order}")

    msa, msa_names = progressive(sequences, names, merge_order)