import numpy as np

# Define the input matrix
input_matrix = np.array([
    [8, 0, 1, 4, 5],
    [1, 0, 5, 0, 2],
    [1, 4, 7, 9, 7],
    [7, 6, 6, 3, 1],
    [0, 8, 6, 0, 0],
    [5, 7, 5, 3, 9],
    [5, 3, 7, 3, 9],
    [3, 8, 5, 0, 6]
])

# Kernel size and stride
kernel_size = 2
stride = 2

# Calculate the dimensions of the output matrix
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the 2x2 submatrix
        submatrix = input_matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Calculate the average of the submatrix
        avg_value = np.mean(submatrix)
        # Assign the average to the output matrix
        output_matrix[i, j] = avg_value

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))