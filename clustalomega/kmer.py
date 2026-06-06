def get_kmers(sqeuence,k = 3):
    """Diziden tüm k-merleri çıkarır"""
    kmers = []
    for i in range(len(sqeuence)-k+1):
        kmers.append(sqeuence[i:i+k]) 
    return kmers 

def kmer_distance(seq1,seq2,k=3):
    """İki dizi arasındaki k-mer mesafesini hesaplar"""
    kmers1=set(get_kmers(seq1,k))
    kmers2=set(get_kmers(seq2,k))

    ortak=len(kmers1 & kmers2) 
    toplam=len(kmers1 | kmers2) 

    if toplam == 0:
        return 0.0 
    
    benzerlik = ortak / toplam 
    return 1- benzerlik 

def distance_matrix(sequences, k=3):
    """Tüm diziler için mesafe matrisi oluşturur"""
    n = len(sequences)
    matrix = []
    
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(0.0)
            else:
                dist = kmer_distance(sequences[i], sequences[j], k)
                row.append(dist)
        matrix.append(row)
    
    return matrix
