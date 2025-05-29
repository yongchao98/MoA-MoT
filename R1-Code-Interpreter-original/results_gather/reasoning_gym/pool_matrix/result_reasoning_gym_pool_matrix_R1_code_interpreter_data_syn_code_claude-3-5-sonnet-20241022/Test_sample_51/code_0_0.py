# Input matrix
matrix = [
    [2, 7, 8, 7, 0, 0, 8, 3, 4],
    [8, 4, 0, 5, 0, 6, 2, 4, 7],
    [5, 0, 9, 9, 5, 9, 0, 4, 6],
    [4, 2, 6, 7, 0, 8, 2, 7, 8],
    [0, 6, 4, 1, 9, 4, 9, 4, 5],
    [1, 5, 3, 7, 9, 3, 9, 5, 2],
    [5, 0, 3, 7, 2, 9, 7, 5, 5],
    [4, 4, 6, 5, 9, 9, 2, 8, 3],
    [4, 7, 2, 4, 8, 1, 9, 7, 9]
]

kernel_size = 1

# Function to perform average pooling
def average_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    result = []
    
    for i in range(0, rows, kernel_size):
        row = []
        for j in range(0, cols, kernel_size):
            # For kernel size 1, just take the element itself
            pool_sum = matrix[i][j]
            pool_avg = pool_sum / (kernel_size * kernel_size)
            row.append(f"{pool_avg:.2f}")
        result.append(row)
    
    # Print in the required format
    print("\n".join([" ".join(row) for row in result]))

average_pooling(matrix, kernel_size)