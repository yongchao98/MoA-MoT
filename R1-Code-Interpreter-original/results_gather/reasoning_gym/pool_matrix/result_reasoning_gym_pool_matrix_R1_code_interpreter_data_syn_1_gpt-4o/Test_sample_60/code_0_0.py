import numpy as np

# Input matrix
matrix = np.array([
    [3, 2, 6, 7, 1, 4, 9, 0, 0, 2],
    [7, 2, 2, 1, 1, 4, 1, 1, 8, 9],
    [4, 9, 9, 6, 1, 7, 1, 4, 5, 2],
    [0, 2, 5, 1, 3, 4, 8, 1, 3, 9],
    [2, 1, 0, 0, 2, 0, 0, 8, 3, 9],
    [0, 5, 0, 4, 5, 5, 1, 0, 7, 3],
    [5, 4, 2, 2, 3, 8, 3, 2, 3, 6],
    [4, 6, 0, 6, 3, 8, 7, 2, 6, 1]
])

# Kernel size and stride
kernel_size = 3
stride = 3

# Calculate output dimensions
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the current submatrix
        submatrix = matrix[i*stride:(i+1)*stride, j*stride:(j+1)*stride]
        # Find the maximum value in the submatrix
        output_matrix[i, j] = np.max(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))