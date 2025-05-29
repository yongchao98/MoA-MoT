import numpy as np

# Define the input matrix
matrix = np.array([
    [2, 7, 3, 6, 0],
    [5, 6, 2, 5, 0],
    [4, 2, 9, 0, 7],
    [2, 7, 9, 5, 3],
    [5, 7, 9, 9, 1]
])

# Define the kernel size
kernel_size = 3

# Calculate the dimensions of the output matrix
output_rows = (matrix.shape[0] - kernel_size) // kernel_size + 1
output_cols = (matrix.shape[1] - kernel_size) // kernel_size + 1

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