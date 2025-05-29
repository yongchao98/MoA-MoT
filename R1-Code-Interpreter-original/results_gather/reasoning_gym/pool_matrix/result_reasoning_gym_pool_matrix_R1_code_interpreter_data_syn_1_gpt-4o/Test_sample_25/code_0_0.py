import numpy as np

# Define the input matrix
input_matrix = np.array([
    [6, 5, 6, 8, 2, 8, 8],
    [1, 7, 9, 5, 8, 8, 1]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = (input_matrix.shape[0] - kernel_size) // stride + 1
output_cols = (input_matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the current pooling region
        region = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Calculate the average of the current region
        output_matrix[i, j] = np.mean(region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))