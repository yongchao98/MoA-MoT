import numpy as np

# Define the input matrix
input_matrix = np.array([
    [5, 2, 2, 6, 7, 0],
    [9, 0, 2, 7, 4, 5]
])

# Kernel size and stride
kernel_size = 1
stride = 1

# Since kernel size is 1, the output is the same as the input
output_matrix = input_matrix

# Print the output matrix
print(output_matrix)