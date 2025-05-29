import numpy as np

# Define the input matrix
matrix = np.array([
    [3, 0, 7, 9, 5],
    [3, 4, 9, 1, 5],
    [0, 6, 4, 4, 2],
    [0, 3, 3, 4, 8],
    [1, 1, 6, 0, 3]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the output matrix dimensions
output_rows = (matrix.shape[0] - kernel_size) // stride + 1
output_cols = (matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 3x3 region
        region = matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average of the region
        output_matrix[i, j] = np.mean(region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))