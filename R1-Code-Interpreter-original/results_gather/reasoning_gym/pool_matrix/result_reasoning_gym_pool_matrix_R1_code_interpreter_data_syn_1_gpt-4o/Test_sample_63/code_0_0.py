import numpy as np

# Define the input matrix
input_matrix = np.array([
    [6, 9, 3, 6, 5, 2, 4, 7, 2, 1],
    [0, 4, 7, 2, 8, 1, 4, 3, 7, 3],
    [4, 3, 4, 2, 5, 4, 1, 8, 4, 6],
    [2, 6, 1, 5, 3, 6, 3, 6, 1, 0],
    [7, 3, 4, 1, 3, 3, 2, 8, 0, 3],
    [7, 3, 1, 7, 9, 4, 8, 4, 4, 5],
    [4, 3, 2, 6, 1, 4, 7, 4, 8, 7],
    [2, 0, 1, 0, 1, 8, 4, 3, 8, 8],
    [3, 0, 3, 4, 0, 3, 7, 4, 4, 8],
    [3, 1, 9, 8, 9, 2, 3, 4, 3, 9]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
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