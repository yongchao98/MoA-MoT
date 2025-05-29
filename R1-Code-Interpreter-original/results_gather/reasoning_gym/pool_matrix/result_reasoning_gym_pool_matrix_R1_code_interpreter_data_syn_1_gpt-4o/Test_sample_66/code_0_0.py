import numpy as np

# Define the input matrix
matrix = np.array([
    [4, 9],
    [1, 4],
    [7, 6],
    [9, 2],
    [7, 2],
    [9, 9],
    [5, 2]
])

# Define the kernel size
kernel_size = 2

# Calculate the number of submatrices
num_submatrices = len(matrix) // kernel_size

# Initialize the output list
output = []

# Perform max pooling
for i in range(num_submatrices):
    submatrix = matrix[i*kernel_size:(i+1)*kernel_size, :]
    max_value = np.max(submatrix)
    output.append(max_value)

# Print the output matrix
print(output)