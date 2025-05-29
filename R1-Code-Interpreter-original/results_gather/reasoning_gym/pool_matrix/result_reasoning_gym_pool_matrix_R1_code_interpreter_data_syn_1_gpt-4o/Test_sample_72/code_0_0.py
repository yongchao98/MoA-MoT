import numpy as np

# Define the input matrix
matrix = np.array([
    [5, 1, 9, 6, 8],
    [8, 4, 0, 3, 7]
])

# Define kernel size and stride
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
        # Extract the current pooling region
        pooling_region = matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average of the pooling region
        output_matrix[i, j] = np.mean(pooling_region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))