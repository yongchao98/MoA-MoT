import numpy as np

# Define the input matrix
matrix = np.array([
    [2, 5, 8, 4, 1],
    [4, 3, 0, 2, 4],
    [6, 3, 9, 4, 9]
])

# Kernel size and stride
kernel_size = 3
stride = 3

# Perform max pooling
max_value = np.max(matrix[:kernel_size, :kernel_size])
print([[max_value]])