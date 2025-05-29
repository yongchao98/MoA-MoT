import numpy as np

# Define the input matrix
input_matrix = np.array([
    [0, 9, 2, 9, 8, 0],
    [8, 5, 9, 2, 2, 3],
    [1, 1, 8, 3, 4, 1],
    [8, 5, 9, 1, 3, 6]
])

# Kernel size
kernel_size = 1

# Since the kernel size is 1, the output matrix is the same as the input matrix
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)