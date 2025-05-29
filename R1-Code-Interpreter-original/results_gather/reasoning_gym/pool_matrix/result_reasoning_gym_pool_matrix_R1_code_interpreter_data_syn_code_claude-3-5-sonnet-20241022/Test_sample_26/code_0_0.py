import numpy as np

# Input matrix
matrix = [
    [3, 0, 7, 9, 5],
    [3, 4, 9, 1, 5],
    [0, 6, 4, 4, 2],
    [0, 3, 3, 4, 8],
    [1, 1, 6, 0, 3]
]

kernel_size = 3

# Convert to numpy array
matrix = np.array(matrix)

# Calculate output dimensions
output_rows = matrix.shape[0] // kernel_size
output_cols = matrix.shape[1] // kernel_size

# Initialize output matrix
output = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the region
        region = matrix[i*kernel_size:(i+1)*kernel_size, 
                       j*kernel_size:(j+1)*kernel_size]
        # Calculate average
        output[i,j] = np.mean(region)

# Format output to 2 decimal places
formatted_output = np.round(output, 2)
print(formatted_output)