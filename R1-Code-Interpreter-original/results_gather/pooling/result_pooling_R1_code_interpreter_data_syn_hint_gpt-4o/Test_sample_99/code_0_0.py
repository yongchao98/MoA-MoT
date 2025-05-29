import numpy as np

# Define the 6x6 matrix
matrix = np.array([
    [3, -10, -3, -1, -10, 8],
    [9, 12, -10, 13, -2, 4],
    [16, 7, 20, 1, 6, 7],
    [6, 11, 11, 3, 19, 1],
    [-2, 1, 11, -4, 8, -4],
    [20, -1, -1, -10, -7, 12]
])

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_dim = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_dim, output_dim))

# Perform average pooling
for i in range(output_dim):
    for j in range(output_dim):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the current window
        window_average = np.mean(window)
        # Store the average in the output matrix
        output_matrix[i, j] = window_average

# Print the resulting matrix
print(output_matrix)