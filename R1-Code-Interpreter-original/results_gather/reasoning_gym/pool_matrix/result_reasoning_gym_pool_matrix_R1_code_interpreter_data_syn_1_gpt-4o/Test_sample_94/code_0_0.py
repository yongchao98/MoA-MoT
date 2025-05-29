import numpy as np

# Define the input matrix
input_matrix = np.array([
    [6, 1, 3, 0, 4, 9, 6, 1, 8],
    [2, 9, 6, 0, 6, 2, 0, 6, 3]
])

# Kernel size
kernel_size = 1

# Since the kernel size is 1, the output matrix is the same as the input matrix
output_matrix = input_matrix

# Print the output matrix with elements formatted to 2 decimal places
formatted_output = np.array2string(output_matrix, formatter={'float_kind':lambda x: "%.2f" % x})
print(formatted_output)