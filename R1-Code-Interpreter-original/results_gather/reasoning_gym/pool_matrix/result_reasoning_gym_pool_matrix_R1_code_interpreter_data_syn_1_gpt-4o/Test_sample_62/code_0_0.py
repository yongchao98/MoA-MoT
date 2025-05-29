import numpy as np

# Define the input matrix
input_matrix = np.array([
    [2, 0, 5, 3],
    [2, 6, 9, 8],
    [3, 6, 0, 8],
    [7, 0, 2, 7],
    [9, 1, 8, 7],
    [0, 8, 8, 8],
    [3, 6, 2, 7],
    [2, 5, 0, 1],
    [0, 1, 7, 3],
    [1, 2, 4, 2]
])

# Kernel size
kernel_size = 1

# Since the kernel size is 1, the output is the same as the input
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)