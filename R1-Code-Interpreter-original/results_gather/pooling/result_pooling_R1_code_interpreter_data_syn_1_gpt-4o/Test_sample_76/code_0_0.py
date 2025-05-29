import numpy as np

# Original 5x5 matrix
matrix = np.array([
    [15, 5, -4, 6, 7],
    [2, -4, -5, 4, -8],
    [-3, 18, -6, -9, -5],
    [-8, -6, -2, -6, 15],
    [4, 1, -2, -7, 13]
])

# Dimensions of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_size, output_size))

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current 3x3 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the result in the output matrix
        output_matrix[i, j] = window_average

# Print the resulting matrix
print(output_matrix)