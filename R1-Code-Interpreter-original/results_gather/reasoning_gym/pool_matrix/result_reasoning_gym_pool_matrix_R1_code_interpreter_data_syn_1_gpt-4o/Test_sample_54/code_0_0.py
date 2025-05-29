import numpy as np

# Define the input matrix
input_matrix = np.array([
    [4, 0, 3, 8, 9, 9, 9, 1, 7, 8],
    [9, 4, 4, 1, 0, 7, 0, 5, 8, 0]
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
        # Compute the average and store it in the output matrix
        output_matrix[i, j] = np.mean(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))