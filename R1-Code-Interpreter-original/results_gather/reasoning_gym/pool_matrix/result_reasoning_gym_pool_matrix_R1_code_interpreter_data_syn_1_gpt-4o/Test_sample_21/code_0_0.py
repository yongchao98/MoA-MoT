import numpy as np

# Define the input matrix
matrix = np.array([
    [9, 5, 5, 6, 4, 0, 8, 9, 5],
    [2, 9, 4, 3, 6, 7, 9, 9, 7],
    [3, 6, 0, 6, 2, 0, 3, 5, 0],
    [5, 5, 5, 0, 7, 9, 6, 6, 4],
    [8, 8, 1, 8, 2, 0, 2, 6, 1]
])

# Define kernel size
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