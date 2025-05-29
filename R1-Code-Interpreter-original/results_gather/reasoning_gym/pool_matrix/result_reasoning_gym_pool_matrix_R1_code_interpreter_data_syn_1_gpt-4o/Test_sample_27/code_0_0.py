import numpy as np

# Define the input matrix
matrix = np.array([
    [1, 7, 5],
    [9, 0, 3],
    [0, 3, 7],
    [9, 0, 9],
    [4, 5, 2]
])

# Kernel size
kernel_size = 3

# Perform max pooling
max_value = np.max(matrix)

# Print the result
print([[max_value]])