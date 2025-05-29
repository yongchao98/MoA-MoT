import numpy as np

# Define the input matrix
matrix = np.array([
    [6, 0, 2, 6, 5],
    [9, 6, 3, 4, 3],
    [3, 0, 4, 1, 6],
    [6, 0, 5, 8, 7],
    [0, 8, 9, 2, 2],
    [5, 8, 3, 9, 5],
    [9, 9, 4, 2, 0],
    [7, 0, 4, 0, 5],
    [0, 2, 0, 2, 7],
    [1, 5, 7, 7, 6]
])

# Define the kernel size
kernel_size = 3

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 3x3 region
        region = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Find the maximum value in the region
        output_matrix[i, j] = np.max(region)

# Print the output matrix
print(output_matrix)