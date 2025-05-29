import numpy as np

# Define the input matrix
input_matrix = np.array([
    [2, 8, 1],
    [9, 3, 5],
    [1, 4, 5],
    [1, 9, 8]
])

# Kernel size and stride
kernel_size = 1
stride = 1

# Output matrix is the same as input matrix for kernel size 1
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)