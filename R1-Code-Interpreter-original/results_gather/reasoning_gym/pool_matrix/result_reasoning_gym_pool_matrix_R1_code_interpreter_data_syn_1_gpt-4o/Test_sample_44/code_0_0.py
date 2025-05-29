import numpy as np

# Define the input matrix
matrix = np.array([
    [1, 7, 7, 0, 1, 2, 8, 7, 2],
    [8, 4, 5, 0, 3, 1, 4, 1, 2],
    [3, 6, 5, 2, 4, 5, 8, 6, 0],
    [4, 6, 4, 2, 3, 4, 8, 8, 8],
    [4, 9, 9, 9, 8, 5, 8, 6, 9],
    [4, 8, 7, 8, 7, 8, 3, 9, 8],
    [7, 2, 5, 5, 9, 0, 4, 4, 2],
])

# Kernel size
k = 2

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // k
output_cols = matrix.shape[1] // k

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the 2x2 submatrix
        submatrix = matrix[i*k:(i+1)*k, j*k:(j+1)*k]
        # Calculate the average and assign it to the output matrix
        output_matrix[i, j] = np.mean(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))