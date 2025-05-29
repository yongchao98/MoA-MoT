# Input matrix
matrix = [
    [1, 7, 2, 9, 0, 4, 4, 5, 5, 6],
    [9, 6, 0, 8, 9, 9, 0, 2, 2, 1],
    [1, 9, 1, 9, 4, 1, 3, 5, 6, 8],
    [7, 5, 5, 0, 8, 4, 9, 3, 9, 2],
    [5, 0, 6, 0, 4, 1, 5, 4, 2, 6],
    [4, 2, 9, 5, 7, 8, 9, 8, 9, 0]
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
        result.append(row)
    
    # Format output
    output = "\n".join([" ".join(row) for row in result])
    print(output)

average_pooling(matrix, 1)