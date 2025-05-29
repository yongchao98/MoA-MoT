import numpy as np

# Define the input matrix
matrix = np.array([
    [0, 3, 9, 0, 0],
    [0, 6, 2, 8, 9],
    [5, 0, 1, 1, 0],
    [8, 7, 9, 7, 1],
    [9, 2, 9, 1, 5],
    [0, 5, 5, 2, 2],
    [3, 6, 4, 0, 1],
    [7, 3, 1, 5, 5]
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