import numpy as np

# Input matrix
matrix = np.array([
    [2, 3, 7],
    [4, 5, 2],
    [0, 3, 0],
    [2, 2, 9],
    [0, 0, 8]
])

def max_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    result_rows = (rows + kernel_size - 1) // kernel_size
    result_cols = (cols + kernel_size - 1) // kernel_size
    result = np.zeros((result_rows, result_cols))
    
    for i in range(0, rows, kernel_size):
        for j in range(0, cols, kernel_size):
            # Define the region for pooling
            region = matrix[i:min(i+kernel_size, rows), j:min(j+kernel_size, cols)]
            result[i//kernel_size, j//kernel_size] = np.max(region)
    
    return result

# Perform max pooling
result = max_pooling(matrix, 2)

# Format output to 2 decimal places
output = np.round(result, 2)
print(output)