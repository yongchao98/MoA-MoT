import numpy as np

# Input matrix
matrix = np.array([
    [2, 2, 2, 6],
    [9, 7, 4, 6],
    [8, 8, 8, 9],
    [4, 1, 2, 8],
    [8, 4, 7, 5],
    [8, 4, 6, 0]
])

kernel_size = 3
rows, cols = matrix.shape
output_rows = (rows + kernel_size - 1) // kernel_size
output_cols = (cols + kernel_size - 1) // kernel_size

# Initialize output matrix
output = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Define the region for pooling
        r_start = i * kernel_size
        r_end = min(r_start + kernel_size, rows)
        c_start = j * kernel_size
        c_end = min(c_start + kernel_size, cols)
        
        # Get maximum value in the region
        region = matrix[r_start:r_end, c_start:c_end]
        output[i, j] = np.max(region)

# Format output
output_str = "\n".join([" ".join([f"{x:.2f}" for x in row]) for row in output])
print(f"<<<\n{output_str}\n>>>")