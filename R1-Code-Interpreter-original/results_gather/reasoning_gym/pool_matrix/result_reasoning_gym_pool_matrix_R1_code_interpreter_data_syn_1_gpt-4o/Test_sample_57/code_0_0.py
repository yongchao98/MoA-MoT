import numpy as np

# Define the input matrix
input_matrix = np.array([
    [3, 4, 6, 5, 8, 2],
    [1, 7, 8, 0, 7, 6],
    [4, 1, 6, 8, 6, 6]
])

# Kernel size and stride
kernel_size = 1
stride = 1

# Output matrix is the same as input matrix for kernel size 1
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)