from .kmer import distance_matrix
from .tree import upgma
from . import progressive

def fasta_oku(dosya):
    """FASTA dosyasından dizileri okur"""
    sequences = []
    names = []
    current_name = None
    current_seq = ""

    with open(dosya, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    names.append(current_name)
                    sequences.append(current_seq)
                current_name = line[1:]
                current_seq = ""
            else:
                current_seq += line

    if current_name:
        names.append(current_name)
        sequences.append(current_seq)

    return sequences, names

def manuel_giris():
    """Kullanıcıdan manuel dizi girişi alır"""
    sequences = []
    names = []

    n = int(input("Kaç dizi gireceksiniz? "))
    for i in range(n):
        name = input(f"\nDizi {i+1} adı: ")
        seq  = input(f"Dizi {i+1} sekansı: ").upper()
        names.append(name)
        sequences.append(seq)

    return sequences, names

def main():
    print("=" * 50)
    print("   CLUSTAL OMEGA — Çoklu Dizi Hizalama")
    print("=" * 50)


    print("\nGiriş yöntemi seçin:")
    print("1. Manuel giriş")
    print("2. FASTA dosyası")
    secim = input("Seçiminiz (1/2): ")

    if secim == "2":
        dosya = input("FASTA dosya yolu: ")
        sequences, names = fasta_oku(dosya)
    else:
        sequences, names = manuel_giris()

  
    k = input("\nK-mer uzunluğu (varsayılan 3): ")
    k = int(k) if k else 3

    print("\n" + "=" * 50)
    print("Girilen diziler:")
    for name, seq in zip(names, sequences):
        print(f"  {name}: {seq}")

    print("\n" + "=" * 50)
    matrix = distance_matrix(sequences, k)

    tree, merge_order = upgma(matrix, names)

    msa, msa_names = progressive(sequences, names, merge_order)

    print("\n" + "=" * 50)
    print("SONUÇ:")
    for name, seq in zip(msa_names, msa):
        print(f"  {name:<15} {seq}")
    print("=" * 50)

if __name__ == "__main__":
    main()