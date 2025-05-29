import numpy as np

# Define the input matrix
matrix = np.array([
    [3, 0, 7, 9, 5],
    [3, 4, 9, 1, 5],
    [0, 6, 4, 4, 2],
    [0, 3, 3, 4, 8],
    [1, 1, 6, 0, 3]
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
        # Extract the current 3x3 block
        block = matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average of the block
        block_average = np.mean(block)
        # Assign the average to the output matrix
        output_matrix[i, j] = block_average

# Format the output matrix to 2 decimal places
output_matrix = np.round(output_matrix, 2)

print(output_matrix)