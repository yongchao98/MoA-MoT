import numpy as np

# Define the input matrix
matrix = np.array([
    [8, 0, 3, 1, 8, 0, 6, 2],
    [1, 0, 7, 7, 6, 3, 9, 6],
    [4, 4, 0, 0, 0, 4, 0, 1],
    [3, 9, 6, 3, 9, 6, 0, 9],
    [1, 8, 0, 7, 8, 1, 6, 1],
    [9, 7, 0, 7, 5, 0, 9, 4],
    [9, 3, 6, 7, 1, 4, 6, 1],
    [3, 4, 4, 8, 4, 5, 2, 5],
    [9, 4, 4, 5, 0, 0, 0, 5]
])

# Define the kernel size
kernel_size = 2

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        submatrix = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        output_matrix[i, j] = np.mean(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))