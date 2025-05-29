import numpy as np

# Input matrix
matrix = np.array([
    [4, 3, 1, 3, 0, 4, 3, 8, 7, 7],
    [6, 9, 3, 7, 3, 3, 6, 5, 4, 5],
    [9, 1, 8, 7, 4, 5, 3, 0, 4, 9],
    [2, 8, 8, 6, 2, 0, 3, 4, 8, 3],
    [2, 2, 1, 2, 2, 9, 8, 1, 8, 9],
    [4, 2, 4, 6, 7, 5, 5, 6, 2, 5],
    [1, 8, 9, 1, 8, 0, 9, 3, 5, 9],
    [5, 0, 8, 0, 4, 2, 9, 7, 6, 6]
])

# Kernel size
kernel_size = 3

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the 3x3 submatrix
        submatrix = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Calculate the average and store it in the output matrix
        output_matrix[i, j] = np.mean(submatrix)

# Round the output matrix to 2 decimal places
output_matrix = np.round(output_matrix, 2)

# Print the output matrix
print(output_matrix)