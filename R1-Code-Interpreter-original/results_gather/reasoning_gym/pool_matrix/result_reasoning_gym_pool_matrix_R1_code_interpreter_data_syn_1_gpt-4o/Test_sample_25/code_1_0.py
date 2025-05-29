import numpy as np

# Define the input matrix
input_matrix = np.array([
    [6, 5, 6, 8, 2, 8, 8],
    [1, 7, 9, 5, 8, 8, 1]
])

# Define kernel size and stride
kernel_size = 3
stride = 3

# Calculate the dimensions of the output matrix
output_cols = input_matrix.shape[1] // kernel_size

# Initialize the output matrix
output_matrix = np.zeros((input_matrix.shape[0], output_cols))

# Perform average pooling along the columns
for i in range(input_matrix.shape[0]):
    for j in range(output_cols):
        # Define the current 1x3 region
        region = input_matrix[i, j*stride:(j+1)*stride]
        # Calculate the average of the region
        region_average = np.mean(region)
        # Assign the average to the output matrix
        output_matrix[i, j] = round(region_average, 2)

# Print the output matrix
print(output_matrix)