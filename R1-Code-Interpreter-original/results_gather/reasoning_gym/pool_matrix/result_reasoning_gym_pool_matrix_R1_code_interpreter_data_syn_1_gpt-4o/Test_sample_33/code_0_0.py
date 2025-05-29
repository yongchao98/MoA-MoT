import numpy as np

# Define the input matrix
matrix = np.array([
    [5, 8, 5, 7, 8, 6],
    [6, 1, 8, 5, 1, 8],
    [9, 1, 7, 6, 8, 3],
    [8, 4, 0, 2, 3, 1],
    [6, 5, 0, 6, 3, 6],
    [1, 7, 7, 2, 4, 5],
    [8, 7, 2, 9, 0, 4],
    [2, 0, 4, 4, 8, 4]
])

# Define kernel size
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