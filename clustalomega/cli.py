from . import align

def main():
    print("=" * 50)
    print("   CLUSTAL OMEGA — Çoklu Dizi Hizalama")
    print("=" * 50)

    print("\nGiriş yöntemi:")
    print("1. Manuel giriş")
    print("2. FASTA dosyası")
    secim = input("Seçiminiz (1/2): ")

    sequences = []
    names = []

    if secim == "2":
        dosya = input("FASTA dosya yolu: ")
        with open(dosya, 'r') as f:
            current_name = None
            current_seq = ""
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
    else:
        n = int(input("Kaç dizi? "))
        for i in range(n):
            name = input(f"Dizi {i+1} adı: ")
            seq  = input(f"Dizi {i+1} sekansı: ").upper()
            names.append(name)
            sequences.append(seq)

    k = input("\nK-mer uzunluğu (varsayılan 3): ")
    k = int(k) if k else 3

    msa, msa_names = align(sequences, names, k)

    print("\n=== SONUÇ ===")
    for name, seq in zip(msa_names, msa):
        print(f"{name:<15} {seq}")

if __name__ == "__main__":
    main()