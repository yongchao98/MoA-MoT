import numpy as np

# Define the input matrix
matrix = np.array([
    [2, 2, 2, 6],
    [9, 7, 4, 6],
    [8, 8, 8, 9],
    [4, 1, 2, 8],
    [8, 4, 7, 5],
    [8, 4, 6, 0]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the current 3x3 region
        region = matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Find the maximum value in the region
        output_matrix[i, j] = np.max(region)

# Print the output matrix
print(output_matrix)