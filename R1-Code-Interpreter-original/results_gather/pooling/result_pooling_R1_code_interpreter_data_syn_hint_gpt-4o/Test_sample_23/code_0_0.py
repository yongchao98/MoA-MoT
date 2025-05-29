import numpy as np

# Original matrix
matrix = np.array([
    [7, 2, 4, 13],
    [18, 13, 10, 13],
    [8, 15, 19, 0],
    [-4, 7, 9, -5]
])

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the resulting matrix
result_rows = matrix.shape[0] - window_size + 1
result_cols = matrix.shape[1] - window_size + 1

# Initialize the resulting matrix
result_matrix = np.zeros((result_rows, result_cols), dtype=int)

# Perform max pooling
for i in range(result_rows):
    for j in range(result_cols):
        # Extract the current 2x2 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        result_matrix[i, j] = np.max(window)

# Print the resulting matrix
print(result_matrix)