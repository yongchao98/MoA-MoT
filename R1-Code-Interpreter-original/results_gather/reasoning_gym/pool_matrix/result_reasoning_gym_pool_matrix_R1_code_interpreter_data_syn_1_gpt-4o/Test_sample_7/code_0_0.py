import numpy as np

# Define the input matrix
input_matrix = np.array([
    [9, 9, 9, 0, 4, 6],
    [3, 0, 4, 8, 1, 4],
    [9, 0, 4, 4, 4, 9],
    [6, 5, 7, 8, 8, 8],
    [7, 1, 9, 7, 6, 9]
])

# Define kernel size and stride
kernel_size = 2
stride = 2

# Calculate output dimensions
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the 2x2 submatrix
        submatrix = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average
        avg_value = np.mean(submatrix)
        # Assign the average to the output matrix
        output_matrix[i, j] = avg_value

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))