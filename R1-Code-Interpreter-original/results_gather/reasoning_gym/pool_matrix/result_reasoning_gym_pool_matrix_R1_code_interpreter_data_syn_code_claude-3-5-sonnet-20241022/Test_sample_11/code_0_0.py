import numpy as np

# Input matrix
matrix = np.array([
    [4, 0, 1, 5, 0, 3],
    [1, 2, 7, 0, 3, 2]
])

kernel_size = 3
rows, cols = matrix.shape
output_rows = rows // kernel_size
output_cols = cols // kernel_size
output = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current kernel region
        region = matrix[
            i*kernel_size:(i+1)*kernel_size,
            j*kernel_size:(j+1)*kernel_size
        ]
        # Calculate average of available elements
        output[i,j] = np.mean(region)

# Format output to 2 decimal places
formatted_output = ' '.join([f'{x:.2f}' for x in output.flatten()])
print(formatted_output)