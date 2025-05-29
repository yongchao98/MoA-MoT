import numpy as np

# Original 5x5 matrix
matrix = np.array([
    [0, 12, 3, 11, -7],
    [7, 10, -2, -7, -3],
    [15, 18, 5, 14, 20],
    [-1, 6, 4, 16, 8],
    [-8, 14, -5, 10, 3]
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
        # Store the average in the output matrix
        output_matrix[i, j] = window_average

# Print the output matrix
print(output_matrix)