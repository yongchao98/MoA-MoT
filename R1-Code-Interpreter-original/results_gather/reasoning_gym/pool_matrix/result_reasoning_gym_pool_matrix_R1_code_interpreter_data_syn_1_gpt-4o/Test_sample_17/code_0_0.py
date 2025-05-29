import numpy as np

# Define the input matrix
input_matrix = np.array([
    [8, 5, 6, 4, 2, 5, 5],
    [1, 2, 4, 0, 2, 3, 6],
    [1, 6, 4, 8, 3, 8, 9]
])

# Define kernel size and stride
kernel_size = 2
stride = 2

# Calculate the dimensions of the output matrix
output_rows = (input_matrix.shape[0] - kernel_size) // stride + 1
output_cols = (input_matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the 2x2 submatrix
        submatrix = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Compute the average and assign it to the output matrix
        output_matrix[i, j] = np.mean(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))