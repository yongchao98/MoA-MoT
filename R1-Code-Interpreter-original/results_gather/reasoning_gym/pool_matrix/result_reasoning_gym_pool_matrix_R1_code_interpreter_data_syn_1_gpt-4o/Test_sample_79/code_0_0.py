import numpy as np

# Define the input matrix
matrix = np.array([
    [4, 4, 4, 6, 3, 9, 2],
    [6, 2, 8, 7, 9, 7, 1],
    [0, 8, 2, 5, 6, 8, 5],
    [6, 5, 8, 6, 2, 6, 4],
    [3, 7, 9, 3, 1, 6, 8]
])

# Define the kernel size
kernel_size = 2

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the 2x2 region
        region = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
        # Calculate the average of the region
        region_average = np.mean(region)
        # Assign the average to the output matrix
        output_matrix[i, j] = region_average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))