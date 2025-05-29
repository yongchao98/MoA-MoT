import numpy as np

# Original matrix
matrix = np.array([
    [2, 3, 7],
    [4, 5, 2],
    [0, 3, 0],
    [2, 2, 9],
    [0, 0, 8]
])

# Kernel size
kernel_size = 2

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the submatrix
        submatrix = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Find the max value in the submatrix
        output_matrix[i, j] = np.max(submatrix)

# Print the output matrix
print(output_matrix)