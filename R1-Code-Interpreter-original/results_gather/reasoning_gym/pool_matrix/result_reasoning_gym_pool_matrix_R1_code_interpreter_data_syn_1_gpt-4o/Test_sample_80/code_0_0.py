import numpy as np

# Define the input matrix
matrix = np.array([
    [4, 8, 1, 4, 5, 7],
    [8, 4, 9, 0, 4, 8],
    [6, 0, 3, 5, 6, 5],
    [9, 8, 3, 1, 5, 2],
    [4, 0, 8, 0, 0, 8],
    [7, 4, 2, 8, 2, 1],
    [3, 4, 2, 3, 7, 0],
    [7, 8, 2, 7, 9, 1],
    [7, 8, 9, 6, 5, 9]
])

# Kernel size and stride
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
        # Extract the current 3x3 region
        region = matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Compute the maximum value in the region
        output_matrix[i, j] = np.max(region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))