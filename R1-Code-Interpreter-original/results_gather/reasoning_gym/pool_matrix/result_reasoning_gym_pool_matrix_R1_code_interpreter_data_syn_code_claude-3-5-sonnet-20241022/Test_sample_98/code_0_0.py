import numpy as np

# Input matrix
matrix = [
    [9, 5, 3, 1, 0, 6, 8, 7],
    [0, 5, 7, 1, 6, 1, 9, 7],
    [2, 1, 4, 8, 2, 2, 1, 2],
    [7, 3, 1, 2, 6, 8, 4, 3],
    [8, 5, 2, 9, 3, 5, 8, 2],
    [7, 7, 2, 0, 5, 1, 0, 8],
    [3, 0, 4, 9, 2, 8, 9, 1],
    [4, 6, 5, 7, 5, 5, 3, 3],
    [1, 1, 2, 5, 5, 9, 7, 2],
    [2, 1, 4, 0, 4, 4, 8, 0]
]

# Convert to numpy array
matrix = np.array(matrix)

# Get dimensions
rows, cols = matrix.shape
kernel_size = 2

# Initialize output matrix
output_rows = rows // kernel_size
output_cols = cols // kernel_size
output = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract current 2x2 block
        block = matrix[i*kernel_size:(i+1)*kernel_size, 
                      j*kernel_size:(j+1)*kernel_size]
        # Get max value from block
        output[i,j] = np.max(block)

# Format output to 2 decimal places
output_formatted = np.array2string(output, formatter={'float_kind':lambda x: "%.2f" % x})
print(output_formatted)