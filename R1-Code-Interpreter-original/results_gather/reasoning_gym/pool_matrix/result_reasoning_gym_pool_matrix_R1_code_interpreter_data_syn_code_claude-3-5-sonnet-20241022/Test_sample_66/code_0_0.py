import numpy as np

# Input matrix
matrix = np.array([
    [4, 9],
    [1, 4],
    [7, 6],
    [9, 2],
    [7, 2],
    [9, 9],
    [5, 2]
])

kernel_size = 2

# Calculate output dimensions
rows = matrix.shape[0]
cols = matrix.shape[1]
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Initialize output matrix
output = np.zeros((out_rows, out_cols))

# Perform max pooling
for i in range(out_rows):
    for j in range(out_cols):
        # Extract the current block
        block = matrix[i*kernel_size:(i+1)*kernel_size, 
                      j*kernel_size:(j+1)*kernel_size]
        # Get maximum value from the block
        output[i,j] = np.max(block)

# Format output to 2 decimal places
formatted_output = np.round(output, 2)
print(formatted_output)