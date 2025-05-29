import numpy as np

# Define the input matrix
input_matrix = np.array([
    [5, 6, 3, 3, 5, 8, 5, 7, 8, 8],
    [3, 2, 8, 1, 8, 3, 8, 4, 1, 3],
    [2, 3, 3, 2, 9, 3, 7, 5, 6, 1]
])

# Kernel size and stride
kernel_size = 1
stride = 1

# Since kernel size is 1, the output matrix is the same as the input matrix
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)