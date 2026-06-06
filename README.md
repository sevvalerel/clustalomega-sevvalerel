# clustalomega_sevvalerel

Clustal Omega çoklu dizi hizalama (Multiple Sequence Alignment) algoritmasının sıfırdan Python implementasyonu. Harici kütüphane gerektirmez, sadece Python standart kütüphanesi kullanılır.

## Kurulum

pip install clustalomega-sevvalerel

Python 3.8 veya üzeri gereklidir.

## Kullanım

### Komut satırı

pip install clustalomega-sevvalerel komutuyla kurduktan sonra terminalde clustalomega yazarak programı başlatabilirsiniz. Program size giriş yöntemi sorar: manuel giriş veya FASTA dosyası. Ardından k-mer uzunluğunu girmeniz istenir ve hizalama otomatik olarak gerçekleştirilir.

### Kütüphane olarak

from clustalomega import align

sequences = ["ACGTAGCT", "ACTAGCT", "TTGCAGCT"]
names = ["Insan", "Fare", "Balik"]

msa, msa_names = align(sequences, names, k=3)

for name, seq in zip(msa_names, msa):
    print(f"{name}: {seq}")

Çıktı:
Insan: ACGT--AGCT
Fare:  AC-T--AGCT
Balik: -T-TGCAGCT

## Algoritma

Algoritma 3 ana adımdan oluşur.

Birinci adımda k-mer yöntemiyle her dizi çifti arasındaki mesafe hesaplanır. K-mer, dizideki uzunluğu k olan alt dizilerin sayılmasıyla elde edilen benzerlik ölçüsüdür. Ortak k-mer sayısı ne kadar yüksekse diziler o kadar benzerdir.

İkinci adımda mesafe matrisinden UPGMA algoritmasıyla guide tree oluşturulur. Bu ağaç hangi dizilerin önce hizalanacağını belirler. En yakın iki dizi önce birleştirilir ve bu işlem tüm diziler tek ağaçta toplanana kadar devam eder.

Üçüncü adımda ağaca göre progressive alignment yapılır. Her adımda hizalanmış dizilerden profil HMM oluşturulur ve Viterbi algoritmasıyla yeni dizi bu profile hizalanır. Viterbi üç durum kullanır: Match (harf eşleşti), Insert (diziye harf eklendi) ve Delete (diziden harf silindi).

## FASTA Dosyası Örneği

>Insan
ACGTAGCTGAC
>Fare
ACTAGCTGAT
>Balik
TTGCAGCTGAC

## Geliştirici

Sevval Erel