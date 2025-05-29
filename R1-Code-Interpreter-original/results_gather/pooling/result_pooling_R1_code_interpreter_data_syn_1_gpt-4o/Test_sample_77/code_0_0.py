import numpy as np

# Define the 5x5 matrix
matrix = np.array([
    [-3, -4, -3, -5, 5],
    [-6, 13, 10, 13, 5],
    [20, 15, -6, 8, 9],
    [8, 9, 4, 4, 10],
    [-1, 7, -2, 1, -7]
])

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_dim = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_dim, output_dim))

# Perform the average pooling
for i in range(output_dim):
    for j in range(output_dim):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the result in the output matrix
        output_matrix[i, j] = window_average

# Print the output matrix
print(output_matrix)