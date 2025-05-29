import numpy as np

# Define the input matrix
input_matrix = np.array([
    [0, 4, 3],
    [4, 3, 5],
    [1, 8, 3],
    [7, 6, 6],
    [0, 7, 9]
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
        # Extract the current 2x2 submatrix
        submatrix = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Compute the average of the submatrix
        avg_value = np.mean(submatrix)
        # Assign the average value to the output matrix
        output_matrix[i, j] = avg_value

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))