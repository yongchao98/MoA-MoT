import numpy as np

# Define the input matrix
matrix = np.array([
    [6, 3],
    [7, 4],
    [6, 9]
])

# Perform max pooling with a kernel size of 3
kernel_size = 3
output = np.max(matrix)

# Print the output matrix
print(np.array([[output]]))