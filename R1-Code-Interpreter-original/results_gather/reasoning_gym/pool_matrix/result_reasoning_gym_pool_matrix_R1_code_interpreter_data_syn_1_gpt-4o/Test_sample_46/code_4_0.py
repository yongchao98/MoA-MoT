import numpy as np

# Define the input matrix
matrix = np.array([
    [5, 6, 7],
    [8, 9, 6],
    [8, 4, 4],
    [3, 3, 8],
    [4, 4, 6]
])

# Define the kernel size
kernel_size = 2

# Calculate the dimensions of the output matrix
output_rows = (matrix.shape[0]) // kernel_size
output_cols = (matrix.shape[1]) // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the submatrix
        submatrix = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Find the maximum value in the submatrix
        output_matrix[i, j] = np.max(submatrix)

# Print the output matrix
print(output_matrix)