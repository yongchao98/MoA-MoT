import numpy as np

# Define the input matrix
input_matrix = np.array([
    [7, 3, 4, 8, 2, 8, 1, 7, 1],
    [7, 2, 9, 4, 2, 6, 4, 3, 6],
    [0, 6, 3, 3, 0, 6, 2, 9, 4],
    [2, 0, 4, 5, 0, 4, 6, 6, 4]
])

# Define kernel size and stride
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
        # Calculate the average and assign it to the output matrix
        output_matrix[i, j] = np.mean(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))