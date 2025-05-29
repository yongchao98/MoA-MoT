import numpy as np

# Define the input matrix
input_matrix = np.array([
    [7, 8, 8, 6, 6, 0],
    [9, 6, 5, 9, 1, 1],
    [0, 9, 8, 1, 1, 6],
    [4, 1, 4, 2, 0, 8],
    [0, 2, 8, 0, 5, 3],
    [0, 3, 9, 8, 6, 6],
    [8, 3, 6, 7, 2, 3],
    [0, 6, 7, 2, 3, 7],
    [5, 5, 7, 3, 2, 1],
    [1, 7, 4, 4, 6, 1]
])

# Kernel size
kernel_size = 1

# Since the kernel size is 1, the output matrix is the same as the input matrix
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)