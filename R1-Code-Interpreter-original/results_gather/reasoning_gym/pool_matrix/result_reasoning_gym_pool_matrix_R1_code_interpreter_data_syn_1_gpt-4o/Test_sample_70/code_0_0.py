import numpy as np

# Define the input matrix
input_matrix = np.array([
    [0, 3, 2, 2, 2, 9, 8],
    [8, 9, 7, 4, 7, 0, 6],
    [2, 7, 3, 0, 0, 0, 6],
    [7, 9, 7, 4, 9, 3, 2]
])

# Kernel size
kernel_size = 1

# Since the kernel size is 1, the output matrix is the same as the input matrix
output_matrix = input_matrix

# Print the output matrix with elements formatted to 2 decimal places
print(np.array2string(output_matrix, formatter={'float_kind':lambda x: "%.2f" % x}))