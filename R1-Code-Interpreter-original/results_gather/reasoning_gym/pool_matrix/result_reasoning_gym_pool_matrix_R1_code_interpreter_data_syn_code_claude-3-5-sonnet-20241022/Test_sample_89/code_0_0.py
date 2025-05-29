import numpy as np

# Input matrix
matrix = np.array([
    [2, 5, 8, 4, 1],
    [4, 3, 0, 2, 4],
    [6, 3, 9, 4, 9]
])

kernel_size = 3

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
        # Extract the current kernel region
        region = matrix[i*kernel_size:(i+1)*kernel_size, 
                       j*kernel_size:(j+1)*kernel_size]
        # Get maximum value in the region
        output[i,j] = np.max(region)

# Format output to 2 decimal places
formatted_output = np.round(output, 2)
print(formatted_output)