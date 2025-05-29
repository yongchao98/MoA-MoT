import numpy as np

# Define the input matrix
matrix = np.array([
    [4, 6, 4, 1, 1, 2, 6],
    [6, 9, 2, 5, 1, 2, 3],
    [1, 7, 1, 4, 8, 4, 5],
    [8, 8, 9, 1, 7, 4, 7],
    [8, 5, 5, 6, 4, 2, 7],
    [9, 5, 3, 4, 2, 9, 8],
    [3, 7, 5, 6, 8, 4, 5]
])

# Define kernel size and stride
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
        # Extract the current 3x3 region
        region = matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Find the maximum value in the region
        output_matrix[i, j] = np.max(region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))