import numpy as np

# Define the original matrix
matrix = np.array([
    [2, 12, 10, 5, 7],
    [19, -5, 17, 20, -1],
    [2, 20, 18, 20, 2],
    [0, -9, -1, -3, -5],
    [10, 5, -6, -5, -2]
])

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the resulting matrix
result_dim = len(matrix) - window_size + 1

# Initialize the resulting matrix
result_matrix = np.zeros((result_dim, result_dim), dtype=int)

# Perform max pooling
for i in range(result_dim):
    for j in range(result_dim):
        # Extract the current 3x3 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        max_value = np.max(window)
        # Store the maximum value in the result matrix
        result_matrix[i, j] = max_value

# Print the resulting matrix
print(result_matrix)