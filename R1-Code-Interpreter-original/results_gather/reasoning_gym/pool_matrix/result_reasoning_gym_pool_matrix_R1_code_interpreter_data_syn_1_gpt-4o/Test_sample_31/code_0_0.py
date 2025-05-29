import numpy as np

# Define the input matrix
matrix = np.array([
    [4, 7, 2, 0, 1, 0, 2, 6, 3],
    [5, 8, 6, 8, 5, 8, 0, 9, 3],
    [9, 7, 5, 9, 4, 6, 4, 8, 5],
    [2, 1, 2, 4, 3, 2, 0, 0, 9],
    [3, 6, 6, 5, 6, 9, 8, 2, 1],
    [6, 4, 3, 0, 7, 0, 2, 0, 9],
    [8, 9, 0, 3, 2, 8, 0, 0, 0],
    [6, 7, 8, 5, 7, 9, 3, 6, 3],
    [8, 8, 8, 1, 6, 1, 1, 4, 9]
])

# Define the kernel size
kernel_size = 2

# Calculate the dimensions of the output matrix
output_rows = (matrix.shape[0] - kernel_size) // kernel_size + 1
output_cols = (matrix.shape[1] - kernel_size) // kernel_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 2x2 submatrix
        submatrix = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Find the maximum value in the submatrix
        output_matrix[i, j] = np.max(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))