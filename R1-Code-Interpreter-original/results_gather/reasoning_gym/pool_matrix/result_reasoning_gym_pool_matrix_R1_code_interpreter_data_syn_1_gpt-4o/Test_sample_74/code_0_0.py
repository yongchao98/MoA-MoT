import numpy as np

# Input matrix
matrix = np.array([7, 9, 4, 4, 3, 3, 5, 2, 5, 0, 3, 1, 5, 2, 6, 9, 4])

# Kernel size and stride
kernel_size = 3
stride = 3

# Perform max pooling
output = [np.max(matrix[i:i+kernel_size]) for i in range(0, len(matrix) - kernel_size + 1, stride)]

# Print the output matrix
print(output)