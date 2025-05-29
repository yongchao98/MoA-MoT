# Input matrix
matrix = [
    [5, 6, 1, 2, 1, 1],
    [8, 7, 2, 8, 7, 3],
    [3, 3, 9, 4, 3, 5],
    [9, 3, 3, 9, 6, 0],
    [3, 7, 3, 1, 1, 5],
    [3, 7, 3, 6, 5, 7],
    [3, 8, 5, 7, 3, 5],
    [9, 3, 7, 6, 8, 4],
    [1, 6, 2, 5, 7, 9],
    [5, 4, 7, 2, 2, 5]
]

def average_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    result = []
    
    for i in range(0, rows, kernel_size):
        row = []
        for j in range(0, cols, kernel_size):
            # For kernel size 1, just take the element itself
            avg = matrix[i][j]
            row.append(f"{avg:.2f}")
        result.append(" ".join(row))
    
    print("<<<")
    print("\n".join(result))
    print(">>>")

average_pooling(matrix, 1)