import numpy as np

# Define the input matrix
input_matrix = np.array([
    [0, 0, 1, 4, 6, 5],
    [6, 6, 5, 2, 7, 4],
    [3, 6, 4, 1, 5, 9],
    [0, 6, 3, 2, 9, 3],
    [3, 3, 2, 0, 3, 2],
    [0, 3, 0, 0, 7, 3],
    [8, 7, 4, 4, 0, 0],
    [3, 3, 1, 9, 4, 5],
    [7, 0, 3, 5, 6, 1],
    [0, 0, 0, 0, 0, 0]
])

# Define kernel size
kernel_size = 2

# Calculate output dimensions
output_height = input_matrix.shape[0] // kernel_size
output_width = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_height, output_width))

# Perform average pooling
for i in range(output_height):
    for j in range(output_width):
        # Extract the 2x2 submatrix
        submatrix = input_matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Calculate the average and store it in the output matrix
        output_matrix[i, j] = np.mean(submatrix)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))