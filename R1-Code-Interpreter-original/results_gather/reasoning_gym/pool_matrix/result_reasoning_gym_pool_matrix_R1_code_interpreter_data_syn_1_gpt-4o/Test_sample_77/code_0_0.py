import numpy as np

# Define the input matrix
matrix = np.array([
    [6, 5, 5, 0, 1],
    [4, 9, 3, 0, 2],
    [5, 4, 5, 6, 0],
    [4, 6, 5, 6, 8],
    [0, 7, 4, 3, 6],
    [1, 5, 4, 6, 5]
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

# Print the output matrix
print(output_matrix)