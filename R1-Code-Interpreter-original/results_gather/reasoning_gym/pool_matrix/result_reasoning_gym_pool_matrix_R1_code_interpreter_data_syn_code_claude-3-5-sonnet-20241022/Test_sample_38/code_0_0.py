import numpy as np

# Input matrix
matrix = np.array([
    [6, 2, 8, 0, 4, 5],
    [5, 5, 4, 7, 8, 0],
    [7, 4, 1, 4, 2, 6]
])

kernel_size = 3
rows, cols = matrix.shape
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Initialize output matrix
output = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current region
        region = matrix[
            i*kernel_size:(i+1)*kernel_size,
            j*kernel_size:(j+1)*kernel_size
        ]
        # Find maximum value in the region
        output[i,j] = np.max(region)

# Format output to 2 decimal places
formatted_output = ' '.join([f'{x:.2f}' for x in output.flatten()])
print(formatted_output)