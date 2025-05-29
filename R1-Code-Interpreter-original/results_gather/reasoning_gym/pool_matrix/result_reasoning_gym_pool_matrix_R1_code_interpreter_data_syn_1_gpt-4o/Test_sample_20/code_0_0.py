import numpy as np

# Define the input matrix
matrix = np.array([
    [5, 7],
    [6, 7],
    [0, 5],
    [3, 3],
    [1, 3],
    [1, 5]
])

# Define the kernel size
kernel_size = 3

# Calculate the number of output rows
output_rows = matrix.shape[0] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, matrix.shape[1]))

# Perform max pooling
for i in range(output_rows):
    submatrix = matrix[i*kernel_size:(i+1)*kernel_size, :]
    output_matrix[i] = np.max(submatrix, axis=0)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))