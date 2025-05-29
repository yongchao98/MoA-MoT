import numpy as np

# Define the input matrix
input_matrix = np.array([
    [3, 7, 1],
    [6, 2, 4],
    [2, 9, 5],
    [6, 6, 9]
])

# Define the kernel size
kernel_size = 3

# Calculate the dimensions of the output matrix
output_rows = (input_matrix.shape[0] - kernel_size) // kernel_size + 1
output_cols = (input_matrix.shape[1] - kernel_size) // kernel_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 3x3 region
        region = input_matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Find the maximum value in the region
        output_matrix[i, j] = np.max(region)

# Print the output matrix
print(output_matrix)