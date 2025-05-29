import numpy as np

# Define the input matrix
input_matrix = np.array([
    [2, 2, 6, 1, 3, 9],
    [6, 1, 0, 1, 9, 0]
])

# Kernel size and stride
kernel_size = 1
stride = 1

# Perform max pooling
output_matrix = input_matrix.copy()

# Print the output matrix
print(output_matrix)