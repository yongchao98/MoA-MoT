import numpy as np

# Input matrix
matrix = np.array([
    [6, 2, 0],
    [8, 3, 1],
    [6, 7, 1],
    [4, 5, 7],
    [7, 4, 9],
    [5, 3, 2]
])

kernel_size = 3
rows, cols = matrix.shape

# Calculate output dimensions
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
        output[i, j] = np.max(region)

# Format output
result = ""
for row in output:
    result += " ".join(f"{x:.2f}" for x in row) + "\n"
print(result.strip())