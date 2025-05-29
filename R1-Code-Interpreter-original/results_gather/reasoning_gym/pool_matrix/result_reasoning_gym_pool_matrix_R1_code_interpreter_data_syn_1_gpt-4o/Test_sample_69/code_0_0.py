import numpy as np

# Define the input matrix
matrix = np.array([
    [1, 6, 7, 9, 8, 4, 8, 5, 0, 5],
    [8, 1, 3, 8, 3, 3, 2, 8, 9, 3],
    [7, 0, 9, 7, 9, 8, 4, 3, 3, 7],
    [4, 8, 7, 6, 9, 4, 2, 7, 7, 7],
    [0, 4, 1, 8, 3, 1, 8, 4, 3, 2]
])

# Define the kernel size
kernel_size = 2

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 2x2 region
        region = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Find the maximum value in the region
        output_matrix[i, j] = np.max(region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))