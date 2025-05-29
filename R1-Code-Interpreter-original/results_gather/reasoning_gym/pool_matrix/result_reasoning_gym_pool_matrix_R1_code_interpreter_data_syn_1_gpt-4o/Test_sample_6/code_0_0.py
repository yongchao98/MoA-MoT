import numpy as np

# Define the input matrix
matrix = np.array([
    [1, 5, 3, 4, 2, 1, 7, 9],
    [0, 1, 4, 3, 1, 7, 2, 1],
    [8, 2, 3, 3, 4, 4, 3, 4],
    [7, 5, 3, 7, 4, 9, 2, 3],
    [5, 9, 7, 1, 4, 5, 7, 4],
    [6, 2, 6, 5, 0, 5, 0, 5],
    [3, 3, 6, 5, 2, 5, 6, 7],
    [1, 1, 7, 4, 9, 3, 0, 3]
])

# Kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = (matrix.shape[0] - kernel_size) // stride + 1
output_cols = (matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the current 3x3 region
        region = matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average of the region
        region_average = np.mean(region)
        # Assign the average to the output matrix
        output_matrix[i, j] = region_average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))