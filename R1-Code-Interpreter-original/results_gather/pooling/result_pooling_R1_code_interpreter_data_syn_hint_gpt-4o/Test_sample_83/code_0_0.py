import numpy as np

# Original 5x5 matrix
matrix = np.array([
    [-5, -2, -1, 9, 11],
    [9, 6, 12, 16, 16],
    [11, 14, 18, 0, 16],
    [20, 14, 15, -7, 3],
    [5, 6, -1, 1, 7]
])

# Size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_size, output_size))

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the result in the output matrix
        output_matrix[i, j] = window_average

# Print the output matrix
print(output_matrix)