import numpy as np

# Define the input matrix
input_matrix = np.array([
    [9, 2, 4],
    [3, 9, 4]
])

# Kernel size
kernel_size = 1

# Perform average pooling
output_matrix = np.zeros_like(input_matrix, dtype=float)

for i in range(0, input_matrix.shape[0], kernel_size):
    for j in range(0, input_matrix.shape[1], kernel_size):
        # Extract the current region
        region = input_matrix[i:i+kernel_size, j:j+kernel_size]
        # Compute the average
        output_matrix[i, j] = np.mean(region)

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))