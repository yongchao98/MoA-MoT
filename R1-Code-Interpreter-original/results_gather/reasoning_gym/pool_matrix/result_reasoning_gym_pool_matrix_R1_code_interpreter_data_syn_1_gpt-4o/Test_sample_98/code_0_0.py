import numpy as np

# Define the input matrix
input_matrix = np.array([
    [9, 5, 3, 1, 0, 6, 8, 7],
    [0, 5, 7, 1, 6, 1, 9, 7],
    [2, 1, 4, 8, 2, 2, 1, 2],
    [7, 3, 1, 2, 6, 8, 4, 3],
    [8, 5, 2, 9, 3, 5, 8, 2],
    [7, 7, 2, 0, 5, 1, 0, 8],
    [3, 0, 4, 9, 2, 8, 9, 1],
    [4, 6, 5, 7, 5, 5, 3, 3],
    [1, 1, 2, 5, 5, 9, 7, 2],
    [2, 1, 4, 0, 4, 4, 8, 0]
])

# Define the kernel size
kernel_size = 2

# Calculate the dimensions of the output matrix
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 2x2 submatrix
        submatrix = input_matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Find the maximum value in the submatrix
        output_matrix[i, j] = np.max(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))