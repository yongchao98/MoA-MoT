import numpy as np

# Define the input matrix
input_matrix = np.array([
    [3, 9, 3, 5, 7, 3, 5, 4],
    [8, 2, 3, 4, 4, 4, 4, 3],
    [0, 8, 3, 3, 8, 7, 9, 9],
    [8, 3, 5, 9, 0, 3, 3, 4],
    [7, 1, 9, 5, 1, 1, 7, 7]
])

# Kernel size and stride
kernel_size = 2
stride = 2

# Calculate the dimensions of the output matrix
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 2x2 block
        block = input_matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Calculate the average of the block
        block_average = np.mean(block)
        # Assign the average to the output matrix
        output_matrix[i, j] = block_average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))