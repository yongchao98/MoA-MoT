import numpy as np

# Define the input matrix
input_matrix = np.array([
    [9, 5, 0, 9, 5, 0, 5, 4, 9, 5],
    [5, 4, 4, 3, 6, 4, 6, 9, 1, 4],
    [7, 0, 8, 8, 6, 4, 5, 9, 2, 3],
    [7, 7, 3, 8, 2, 5, 7, 1, 3, 5],
    [8, 9, 6, 5, 7, 5, 1, 4, 1, 1]
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
        # Extract the current 3x3 block
        block = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average of the block
        block_average = np.mean(block)
        # Assign the average to the output matrix
        output_matrix[i, j] = block_average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))