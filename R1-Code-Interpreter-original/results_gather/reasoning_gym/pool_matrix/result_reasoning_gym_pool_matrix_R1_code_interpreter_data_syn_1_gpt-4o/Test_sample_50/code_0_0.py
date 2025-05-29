import numpy as np

# Define the input matrix
input_matrix = np.array([
    [3, 2, 3, 6, 2, 0, 1, 1, 3],
    [3, 9, 5, 4, 0, 1, 9, 7, 1],
    [3, 4, 8, 8, 9, 3, 4, 0, 0],
    [7, 3, 1, 3, 0, 4, 6, 6, 2]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the current 3x3 region
        region = input_matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Find the maximum value in the current region
        max_value = np.max(region)
        # Assign the maximum value to the output matrix
        output_matrix[i, j] = max_value

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))