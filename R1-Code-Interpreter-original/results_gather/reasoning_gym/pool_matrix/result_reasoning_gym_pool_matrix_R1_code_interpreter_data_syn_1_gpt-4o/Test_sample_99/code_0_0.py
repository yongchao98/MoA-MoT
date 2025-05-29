import numpy as np

# Define the input matrix
input_matrix = np.array([
    [3, 7, 9, 9, 5, 0, 2],
    [8, 1, 8, 0, 5, 4, 0],
    [0, 1, 1, 2, 7, 5, 8],
    [7, 5, 0, 7, 0, 3, 8],
    [8, 8, 6, 8, 3, 2, 4],
    [5, 3, 5, 3, 8, 0, 9],
    [1, 8, 2, 2, 6, 5, 6]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_rows = input_matrix.shape[0] // kernel_size
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 3x3 region
        region = input_matrix[i*stride:(i*stride)+kernel_size, j*stride:(j*stride)+kernel_size]
        # Calculate the average of the region
        region_average = np.mean(region)
        # Assign the average to the output matrix
        output_matrix[i, j] = region_average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))