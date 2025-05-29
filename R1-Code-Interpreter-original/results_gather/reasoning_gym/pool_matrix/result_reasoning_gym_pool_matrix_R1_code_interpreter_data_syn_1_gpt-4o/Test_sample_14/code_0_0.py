import numpy as np

# Define the input matrix
input_matrix = np.array([
    [5, 8, 4, 3, 2, 0, 2, 8],
    [0, 7, 6, 2, 2, 2, 2, 6],
    [7, 8, 8, 6, 0, 2, 4, 8]
])

# Kernel size
kernel_size = 1

# Perform average pooling
output_matrix = input_matrix.copy()

# Print the output matrix
print(output_matrix)