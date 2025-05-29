import numpy as np

# Define the input matrix
matrix = np.array([
    [6, 3, 2, 3, 6, 7, 7, 3],
    [1, 5, 2, 4, 0, 9, 9, 4],
    [3, 2, 3, 9, 8, 7, 5, 0],
    [8, 3, 7, 3, 2, 0, 2, 2],
    [9, 8, 6, 8, 1, 8, 1, 6],
    [1, 6, 7, 9, 0, 9, 4, 8]
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
        # Extract the current 3x3 region
        region = matrix[i*stride:(i+1)*stride, j*stride:(j+1)*stride]
        # Find the maximum value in the region
        output_matrix[i, j] = np.max(region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))