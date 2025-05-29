import numpy as np

# Define the input matrix
input_matrix = np.array([
    [7, 8, 8, 2, 1, 0, 6, 1, 3, 9],
    [0, 8, 0, 8, 1, 0, 5, 0, 2, 1],
    [6, 2, 6, 7, 2, 1, 6, 6, 3, 0],
    [0, 1, 0, 8, 0, 7, 8, 0, 5, 8],
    [2, 3, 4, 9, 6, 0, 3, 2, 3, 0],
    [2, 0, 0, 5, 7, 3, 3, 3, 3, 7],
    [6, 8, 9, 1, 6, 5, 4, 8, 3, 7]
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