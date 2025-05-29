import numpy as np

# Define the input matrix
input_matrix = np.array([
    [3, 5, 5],
    [7, 9, 9],
    [3, 6, 9]
])

# Define kernel size and stride
kernel_size = 2
stride = 2

# Calculate the dimensions of the output matrix
output_rows = (input_matrix.shape[0] - kernel_size) // stride + 1
output_cols = (input_matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 2x2 submatrix
        submatrix = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Find the maximum value in the submatrix
        output_matrix[i, j] = np.max(submatrix)

# Print the output matrix
print(output_matrix)