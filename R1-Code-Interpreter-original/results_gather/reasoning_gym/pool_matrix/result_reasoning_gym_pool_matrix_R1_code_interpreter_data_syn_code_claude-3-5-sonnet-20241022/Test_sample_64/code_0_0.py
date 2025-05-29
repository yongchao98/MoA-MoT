import numpy as np

# Input matrix
matrix = [
    [4, 6, 4, 1, 1, 2, 6],
    [6, 9, 2, 5, 1, 2, 3],
    [1, 7, 1, 4, 8, 4, 5],
    [8, 8, 9, 1, 7, 4, 7],
    [8, 5, 5, 6, 4, 2, 7],
    [9, 5, 3, 4, 2, 9, 8],
    [3, 7, 5, 6, 8, 4, 5]
]

kernel_size = 3

def max_pooling(matrix, kernel_size):
    n = len(matrix)
    output_size = (n + kernel_size - 1) // kernel_size
    result = []
    
    for i in range(0, n, kernel_size):
        row = []
        for j in range(0, n, kernel_size):
            # Extract the kernel region
            region = []
            for ki in range(kernel_size):
                for kj in range(kernel_size):
                    if i+ki < n and j+kj < n:
                        region.append(matrix[i+ki][j+kj])
            row.append(max(region))
        result.append(row)
    
    return result

# Perform max pooling
result = max_pooling(matrix, kernel_size)

# Print result in a formatted way
for row in result:
    print(' '.join(f'{x:>2}' for x in row))