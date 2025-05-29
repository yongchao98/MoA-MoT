import numpy as np

# Input matrix
matrix = [
    [4, 7, 2, 0, 1, 0, 2, 6, 3],
    [5, 8, 6, 8, 5, 8, 0, 9, 3],
    [9, 7, 5, 9, 4, 6, 4, 8, 5],
    [2, 1, 2, 4, 3, 2, 0, 0, 9],
    [3, 6, 6, 5, 6, 9, 8, 2, 1],
    [6, 4, 3, 0, 7, 0, 2, 0, 9],
    [8, 9, 0, 3, 2, 8, 0, 0, 0],
    [6, 7, 8, 5, 7, 9, 3, 6, 3],
    [8, 8, 8, 1, 6, 1, 1, 4, 9]
]

# Convert to numpy array
matrix = np.array(matrix)
kernel_size = 2
output_size = matrix.shape[0] // kernel_size

# Initialize output matrix
output = np.zeros((output_size, output_size))

# Perform max pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the 2x2 region
        region = matrix[i*kernel_size:(i+1)*kernel_size, 
                       j*kernel_size:(j+1)*kernel_size]
        # Get maximum value
        output[i,j] = np.max(region)

# Format output to 2 decimal places
formatted_output = np.array2string(output, formatter={'float_kind':lambda x: "%.2f" % x})
print(formatted_output)