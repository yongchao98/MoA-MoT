import numpy as np

# Define the input matrix
matrix = np.array([
    [6, 2, 0],
    [8, 3, 1],
    [6, 7, 1],
    [4, 5, 7],
    [7, 4, 9],
    [5, 3, 2]
])

# Kernel size
kernel_size = 3

# Calculate the number of regions
num_regions = matrix.shape[0] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((num_regions, 1))

# Perform max pooling
for i in range(num_regions):
    region = matrix[i*kernel_size:(i+1)*kernel_size, :]
    output_matrix[i, 0] = np.max(region)

# Print the output matrix
print(output_matrix)