import numpy as np

# Define the input matrix
matrix = np.array([
    [4, 0, 1, 5, 0, 3],
    [1, 2, 7, 0, 3, 2]
])

# Kernel size and stride
kernel_size = 3
stride = 3

# Calculate the number of output columns
output_cols = matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((1, output_cols))

# Perform average pooling
for i in range(output_cols):
    # Extract the current pooling region
    pooling_region = matrix[:, i*stride:(i*stride)+kernel_size]
    # Calculate the average
    average = np.mean(pooling_region)
    # Assign to the output matrix
    output_matrix[0, i] = average

# Print the output matrix with elements rounded to 2 decimal places
print(np.round(output_matrix, 2))