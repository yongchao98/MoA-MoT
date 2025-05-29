import numpy as np

# Define the input matrix
input_matrix = np.array([
    [1, 6, 3, 8, 9, 1],
    [2, 1, 7, 8, 5, 2],
    [8, 6, 3, 2, 5, 0],
    [3, 4, 0, 6, 6, 9],
    [5, 7, 0, 2, 5, 8],
    [7, 9, 4, 0, 7, 0]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 3x3 block
        block = input_matrix[i*stride:(i+1)*stride, j*stride:(j+1)*stride]
        # Calculate the average of the block
        block_average = np.mean(block)
        # Assign the average to the output matrix
        output_matrix[i, j] = block_average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))