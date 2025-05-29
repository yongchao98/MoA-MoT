import numpy as np

# Original 4x4 matrix
matrix = np.array([
    [17, 8, 5, 8],
    [10, -6, -6, 4],
    [0, 4, 16, -9],
    [17, 7, 6, -8]
])

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_size, output_size))

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current 2x2 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the average in the output matrix
        output_matrix[i, j] = window_average

# Print the resulting matrix
print(output_matrix)