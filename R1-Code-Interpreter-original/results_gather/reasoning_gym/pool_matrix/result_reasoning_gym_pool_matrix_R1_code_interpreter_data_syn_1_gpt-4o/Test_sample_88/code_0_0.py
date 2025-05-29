import numpy as np

# Define the input matrix
input_matrix = np.array([
    [7, 0, 8, 4, 1],
    [7, 1, 3, 8, 9],
    [4, 9, 9, 6, 7],
    [9, 1, 9, 4, 5],
    [7, 2, 2, 2, 0],
    [0, 7, 4, 2, 2],
    [4, 8, 5, 6, 4],
    [1, 9, 2, 0, 5],
    [7, 5, 7, 8, 4],
    [2, 1, 6, 9, 5]
])

# Kernel size
kernel_size = 1

# Since the kernel size is 1, the output matrix is the same as the input matrix
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)