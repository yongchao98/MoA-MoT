import numpy as np

# Define the input matrix
matrix = np.array([
    [4, 6, 7],
    [6, 6, 1],
    [4, 3, 5],
    [0, 8, 6],
    [8, 2, 0],
    [1, 4, 2],
    [6, 4, 5],
    [0, 5, 9],
    [6, 2, 6]
])

# Define the kernel size
kernel_size = 3

# Calculate the number of regions
num_regions = matrix.shape[0] // kernel_size

# Initialize the output matrix
output_matrix = []

# Perform max pooling
for i in range(num_regions):
    region = matrix[i*kernel_size:(i+1)*kernel_size, :]
    max_value = np.max(region)
    output_matrix.append(max_value)

# Print the output matrix
print(output_matrix)