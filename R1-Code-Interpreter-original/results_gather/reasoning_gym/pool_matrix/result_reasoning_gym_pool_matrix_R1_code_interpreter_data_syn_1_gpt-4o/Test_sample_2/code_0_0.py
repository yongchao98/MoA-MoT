import numpy as np

# Define the input matrix
matrix = np.array([
    [1, 9, 6, 5, 4, 3, 5],
    [2, 9, 8, 6, 5, 6, 4],
    [5, 0, 5, 4, 3, 3, 5],
    [5, 8, 8, 4, 8, 4, 2],
    [5, 2, 2, 9, 6, 8, 2],
    [6, 0, 5, 4, 0, 4, 4],
    [9, 7, 4, 9, 7, 3, 6],
    [7, 9, 6, 8, 8, 9, 4],
    [4, 0, 1, 1, 5, 9, 3],
    [5, 0, 0, 2, 3, 4, 1]
])

# Kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the submatrix
        submatrix = matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Calculate the average and store it in the output matrix
        output_matrix[i, j] = np.mean(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))