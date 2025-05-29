import numpy as np

# Input matrix
matrix = np.array([
    [3, 1, 8, 3, 1, 1, 1, 3, 9],
    [4, 3, 8, 4, 1, 0, 5, 4, 1],
    [2, 4, 1, 3, 3, 1, 0, 9, 1],
    [6, 3, 6, 9, 1, 5, 2, 6, 9],
    [1, 0, 0, 4, 5, 5, 6, 1, 3],
    [2, 9, 3, 8, 4, 4, 6, 9, 9],
    [5, 2, 6, 0, 1, 4, 4, 4, 2],
    [9, 2, 7, 8, 1, 0, 1, 2, 6],
    [5, 1, 9, 7, 0, 0, 5, 8, 6],
    [4, 8, 1, 1, 0, 1, 7, 6, 3]
])

# Kernel size
k = 2

# Output matrix dimensions
output_rows = matrix.shape[0] // k
output_cols = matrix.shape[1] // k

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the 2x2 submatrix
        submatrix = matrix[i*k:(i+1)*k, j*k:(j+1)*k]
        # Find the maximum value in the submatrix
        output_matrix[i, j] = np.max(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))