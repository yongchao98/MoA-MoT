import numpy as np

# Original matrix
matrix = np.array([
    [1, 9],
    [8, 5],
    [0, 1],
    [4, 1],
    [9, 5],
    [9, 8]
])

# Kernel size
kernel_size = 2

# Perform average pooling
output_matrix = []
for i in range(0, matrix.shape[0], kernel_size):
    if i + kernel_size <= matrix.shape[0]:
        region = matrix[i:i+kernel_size, :]
        avg_values = np.mean(region, axis=0)
        output_matrix.append(avg_values)

# Convert to numpy array and format to 2 decimal places
output_matrix = np.array(output_matrix)
output_matrix = np.round(output_matrix, 2)

print(output_matrix)