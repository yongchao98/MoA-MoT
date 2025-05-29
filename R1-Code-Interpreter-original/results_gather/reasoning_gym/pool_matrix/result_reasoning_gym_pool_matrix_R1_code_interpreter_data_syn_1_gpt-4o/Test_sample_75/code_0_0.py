import numpy as np

# Define the input matrix
matrix = np.array([
    [0, 9, 6, 2, 9, 6, 0, 2],
    [2, 5, 9, 4, 0, 9, 8, 4],
    [3, 8, 0, 5, 6, 9, 7, 3],
    [6, 6, 9, 0, 1, 9, 4, 4],
    [5, 8, 4, 2, 9, 9, 8, 1]
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