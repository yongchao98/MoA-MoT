import numpy as np

# Define the input matrix
input_matrix = np.array([
    [7, 5, 6, 8, 3],
    [7, 9, 3, 5, 9],
    [4, 1, 3, 1, 2]
])

# Kernel size and stride
kernel_size = 1
stride = 1

# Output matrix dimensions
output_height = (input_matrix.shape[0] - kernel_size) // stride + 1
output_width = (input_matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_height, output_width))

# Perform max pooling
for i in range(output_height):
    for j in range(output_width):
        # Extract the current region
        current_region = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Apply max pooling
        output_matrix[i, j] = np.max(current_region)

# Print the output matrix
print(output_matrix)