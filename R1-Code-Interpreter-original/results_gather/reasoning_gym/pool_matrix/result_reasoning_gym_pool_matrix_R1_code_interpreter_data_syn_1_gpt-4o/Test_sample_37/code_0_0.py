import numpy as np

# Define the input matrix
input_matrix = np.array([
    [6, 9, 7],
    [4, 1, 6]
])

# Kernel size and stride
kernel_size = 1
stride = 1

# Calculate the dimensions of the output matrix
output_height = (input_matrix.shape[0] - kernel_size) // stride + 1
output_width = (input_matrix.shape[1] - kernel_size) // stride + 1

# Initialize the output matrix
output_matrix = np.zeros((output_height, output_width))

# Perform average pooling
for i in range(0, output_height):
    for j in range(0, output_width):
        # Extract the current region
        region = input_matrix[i*stride:i*stride+kernel_size, j*stride:j*stride+kernel_size]
        # Compute the average
        output_matrix[i, j] = np.mean(region)

# Print the output matrix
print(np.round(output_matrix, 2))