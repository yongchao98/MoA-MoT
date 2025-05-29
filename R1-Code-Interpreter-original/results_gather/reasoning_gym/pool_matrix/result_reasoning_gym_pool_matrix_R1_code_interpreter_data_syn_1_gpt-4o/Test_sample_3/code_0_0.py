import numpy as np

# Define the input matrix
matrix = np.array([
    [9, 7, 1, 0, 9, 8, 8],
    [8, 9, 7, 7, 8, 0, 9],
    [2, 5, 3, 3, 3, 9, 8],
    [7, 5, 8, 6, 0, 1, 9],
    [9, 0, 2, 2, 5, 8, 3],
    [0, 4, 1, 3, 6, 4, 2],
    [3, 1, 3, 4, 2, 7, 2],
    [4, 7, 4, 7, 9, 4, 7],
    [4, 0, 8, 7, 6, 3, 6],
    [0, 5, 5, 4, 2, 0, 4]
])

# Kernel size
kernel_size = 1

# Perform max pooling
output_matrix = matrix.copy()

# Print the output matrix
print(output_matrix)