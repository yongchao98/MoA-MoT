import numpy as np

# Define the input matrix
input_matrix = np.array([
    [7, 4, 4, 5, 8, 0, 9],
    [7, 5, 3, 0, 6, 4, 6],
    [4, 5, 7, 3, 1, 1, 9],
    [6, 4, 5, 8, 1, 8, 1],
    [6, 3, 4, 3, 9, 0, 8],
    [9, 5, 3, 2, 4, 4, 5]
])

# Kernel size
kernel_size = 1

# Perform average pooling
output_matrix = input_matrix.copy()

# Print the output matrix
print(output_matrix)