import numpy as np

# Define the input matrix
matrix = np.array([
    [5, 2, 7, 1, 0, 8, 7, 1, 9],
    [6, 9, 7, 4, 4, 4, 6, 4, 0],
    [6, 8, 9, 7, 3, 2, 1, 2, 6],
    [9, 8, 1, 7, 7, 2, 7, 7, 8]
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