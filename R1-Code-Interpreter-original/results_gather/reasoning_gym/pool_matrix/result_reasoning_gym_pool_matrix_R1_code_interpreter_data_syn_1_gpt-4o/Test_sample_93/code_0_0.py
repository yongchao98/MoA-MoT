import numpy as np

# Define the input matrix
matrix = np.array([
    [2, 7, 3, 6, 0],
    [5, 6, 2, 5, 0],
    [4, 2, 9, 0, 7],
    [2, 7, 9, 5, 3],
    [5, 7, 9, 9, 1]
])

# Define the kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = (matrix.shape[0] - kernel_size) // stride + 1
output_cols = (matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the current region
        region = matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Find the maximum value in the current region
        output_matrix[i, j] = np.max(region)

# Print the output matrix
print(output_matrix)