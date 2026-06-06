def find_min_distance(matrix, labels):
    """Matristeki en küçük mesafeyi ve indekslerini bulur"""
    min_dist = float('inf')
    min_i = 0
    min_j = 1

    for i in range(len(matrix)):
        for j in range(i+1, len(matrix)):
            if matrix[i][j] < min_dist:
                min_dist = matrix[i][j]
                min_i = i
                min_j = j

    return min_i, min_j, min_dist


def merge_clusters(matrix, labels, i, j):
    """İki kümeyi birleştirir, yeni mesafeleri hesaplar"""
    n = len(matrix)
    new_label = f"({labels[i]},{labels[j]})"
    new_row = []

    for k in range(n):
        if k != i and k != j:
            # Ortalama mesafe hesapla
            new_dist = (matrix[i][k] + matrix[j][k]) / 2
            new_row.append(new_dist)


    new_matrix = []
    new_labels = []

    for k in range(n):
        if k != i and k != j:
            row = []
            for l in range(n):
                if l != i and l != j:
                    row.append(matrix[k][l])
            new_matrix.append(row)
            new_labels.append(labels[k])

    for k in range(len(new_matrix)):
        new_matrix[k].append(new_row[k])
        new_row[k] = new_row[k]

    new_row.append(0.0)
    new_matrix.append(new_row)
    new_labels.append(new_label)

    return new_matrix, new_labels


def upgma(matrix, labels):
    """UPGMA algoritması ile guide tree oluşturur"""
    matrix = [row[:] for row in matrix]
    labels = labels[:]
    merge_order = []  
    print("\n=== UPGMA Guide Tree ===")

    while len(matrix) > 1:
        i, j, dist = find_min_distance(matrix, labels)
        print(f"Birleştiriliyor: {labels[i]} + {labels[j]} (mesafe: {dist:.2f})")
        
        merge_order.append((labels[i], labels[j]))  
        matrix, labels = merge_clusters(matrix, labels, i, j)

    print(f"Guide Tree: {labels[0]}")
    return labels[0], merge_order  


if __name__ == "__main__":
    matrix = [
        [0.0, 0.8, 1.0],
        [0.8, 0.0, 0.8],
        [1.0, 0.8, 0.0]
    ]
    labels = ["Dizi1", "Dizi2", "Dizi3"]

    tree, merge_order = upgma(matrix, labels)
    print(f"Merge sırası: {merge_order}")