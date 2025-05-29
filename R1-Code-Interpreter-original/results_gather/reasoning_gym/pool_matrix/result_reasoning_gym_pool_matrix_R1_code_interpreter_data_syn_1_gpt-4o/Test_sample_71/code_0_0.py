import numpy as np

# Define the input matrix
input_matrix = np.array([
    [7, 9, 3, 9, 6],
    [1, 5, 7, 7, 4],
    [6, 1, 4, 8, 1],
    [4, 5, 0, 1, 6],
    [0, 3, 5, 3, 8],
    [7, 6, 7, 1, 7],
    [8, 3, 5, 1, 8],
    [0, 1, 4, 4, 8],
    [9, 5, 5, 4, 5]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate output dimensions
output_rows = (input_matrix.shape[0] - kernel_size) // stride + 1
output_cols = (input_matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 3x3 region
        region = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average of the region
        region_average = np.mean(region)
        # Assign the average to the output matrix
        output_matrix[i, j] = region_average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))