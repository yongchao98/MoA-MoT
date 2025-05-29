import numpy as np

# Define the input matrix
matrix = np.array([
    [5, 4],
    [0, 2],
    [9, 7],
    [6, 4],
    [9, 7],
    [1, 8],
    [2, 0],
    [5, 6],
    [1, 9],
    [5, 5]
])

# Define the kernel size
kernel_size = 2

# Calculate the number of blocks
num_blocks = matrix.shape[0] // kernel_size

# Initialize the output matrix
output_matrix = []

# Perform max pooling
for i in range(num_blocks):
    block = matrix[i*kernel_size:(i+1)*kernel_size, :]
    max_value = np.max(block)
    output_matrix.append(max_value)

# Convert the output to a numpy array
output_matrix = np.array(output_matrix)

# Print the output matrix
print(output_matrix)