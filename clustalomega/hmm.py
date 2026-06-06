import math

TRANSITION = {
    'MM': 0.7, 'MI': 0.15, 'MD': 0.15,
    'IM': 0.9, 'II': 0.1,
    'DM': 0.9, 'DD': 0.1,
}

PSEUDO = 0.1 
ALFABE = ['A', 'T', 'G', 'C']

def build_profile(aligned_seqs):
    """Hizalanmış dizilerden profil oluşturur"""
    if not aligned_seqs:
        return []
    profile = []
    length = len(aligned_seqs[0])
    for j in range(length):
        counts = {harf: 0 for harf in ALFABE}
        for seq in aligned_seqs:
            if j < len(seq) and seq[j] in counts:
                counts[seq[j]] += 1
        total = len(aligned_seqs) + len(ALFABE) * PSEUDO
        freqs = {harf: (counts[harf] + PSEUDO) / total for harf in ALFABE}
        profile.append(freqs)
    return profile

def emission(profile_col, base):
    """Profil sütununda o harfin frekansını döndürür"""
    if base == '-':
        return 0.0
    return profile_col.get(base, PSEUDO)

def viterbi(profile, sequence):
    n = len(sequence)
    m = len(profile)

    M = [[0.0]*(m+1) for _ in range(n+1)]
    I = [[0.0]*(m+1) for _ in range(n+1)]
    D = [[0.0]*(m+1) for _ in range(n+1)]

    M_back = [['']*(m+1) for _ in range(n+1)]
    I_back = [['']*(m+1) for _ in range(n+1)]
    D_back = [['']*(m+1) for _ in range(n+1)]

    M[0][0] = 1.0

    for j in range(1, m+1):
        v1 = M[0][j-1] * TRANSITION['MD']
        v2 = D[0][j-1] * TRANSITION['DD']
        D[0][j] = max(v1, v2)
        D_back[0][j] = 'M' if v1 >= v2 else 'D'

    for i in range(1, n+1):
        v1 = M[i-1][0] * TRANSITION['MI']
        v2 = I[i-1][0] * TRANSITION['II']
        I[i][0] = max(v1, v2)
        I_back[i][0] = 'M' if v1 >= v2 else 'I'

    for i in range(1, n+1):
        for j in range(1, m+1):
            base = sequence[i-1]
            em = emission(profile[j-1], base)

            m_scores = [
                (M[i-1][j-1] * TRANSITION['MM'], 'M'),
                (I[i-1][j-1] * TRANSITION['IM'], 'I'),
                (D[i-1][j-1] * TRANSITION['DM'], 'D'),
            ]
            best_m = max(m_scores, key=lambda x: x[0])
            M[i][j] = best_m[0] * em
            M_back[i][j] = best_m[1]

            i_scores = [
                (M[i-1][j] * TRANSITION['MI'], 'M'),
                (I[i-1][j] * TRANSITION['II'], 'I'),
            ]
            best_i = max(i_scores, key=lambda x: x[0])
            I[i][j] = best_i[0]
            I_back[i][j] = best_i[1]

            d_scores = [
                (M[i][j-1] * TRANSITION['MD'], 'M'),
                (D[i][j-1] * TRANSITION['DD'], 'D'),
            ]
            best_d = max(d_scores, key=lambda x: x[0])
            D[i][j] = best_d[0]
            D_back[i][j] = best_d[1]

    seq_aligned  = []
    gaps_in_prof = []
    i, j = n, m

    scores = [(M[n][m],'M'), (I[n][m],'I'), (D[n][m],'D')]
    current = max(scores, key=lambda x: x[0])[1]

    while i > 0 or j > 0:
        if i == 0:
            seq_aligned.append('-')
            gaps_in_prof.append(False)
            j -= 1
        elif j == 0:
            seq_aligned.append(sequence[i-1])
            gaps_in_prof.append(True)
            i -= 1
        elif current == 'M':
            seq_aligned.append(sequence[i-1])
            gaps_in_prof.append(False)
            current = M_back[i][j]
            i -= 1
            j -= 1
        elif current == 'I':
            seq_aligned.append(sequence[i-1])
            gaps_in_prof.append(True)
            current = I_back[i][j]
            i -= 1
        else:  # D
            seq_aligned.append('-')
            gaps_in_prof.append(False)
            current = D_back[i][j]
            j -= 1

    seq_aligned  = ''.join(reversed(seq_aligned))
    gaps_in_prof = list(reversed(gaps_in_prof))

    return seq_aligned, gaps_in_prof
