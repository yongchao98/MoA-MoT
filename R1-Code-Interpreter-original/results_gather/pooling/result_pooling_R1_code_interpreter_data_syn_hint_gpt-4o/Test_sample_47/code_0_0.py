import numpy as np

# Original 5x5 matrix
matrix = np.array([
    [17, -2, 14, -1, 19],
    [11, 11, -1, 5, -2],
    [15, 19, 17, 3, -5],
    [-8, 9, 19, -6, 10],
    [12, 14, 11, 18, -10]
])

# Dimensions of the sliding window
window_size = 3

# Calculate the dimensions of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = np.zeros((result_size, result_size), dtype=int)

# Perform max pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        result_matrix[i, j] = np.max(window)

# Print the resulting matrix
print(result_matrix)