import numpy as np

# Define the input matrix
input_matrix = np.array([
    [2, 3, 5, 2, 7, 2, 4, 9],
    [5, 1, 0, 0, 7, 5, 9, 6],
    [5, 0, 6, 3, 1, 3, 6, 0],
    [7, 4, 7, 8, 8, 6, 7, 0],
    [4, 3, 7, 8, 4, 7, 5, 8]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate output dimensions
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the current 3x3 region
        region = input_matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Find the maximum value in the region
        max_value = np.max(region)
        # Assign the max value to the output matrix
        output_matrix[i, j] = max_value

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))