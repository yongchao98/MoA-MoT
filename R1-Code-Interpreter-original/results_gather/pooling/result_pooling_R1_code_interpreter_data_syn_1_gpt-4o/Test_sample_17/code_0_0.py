import numpy as np

# Original matrix
matrix = np.array([
    [9, 3, -4],
    [-6, 4, 12],
    [-10, 12, -4]
])

# Dimensions of the sliding window
window_size = 2

# Resulting matrix dimensions
result_rows = matrix.shape[0] - window_size + 1
result_cols = matrix.shape[1] - window_size + 1

# Initialize the result matrix
result = np.zeros((result_rows, result_cols), dtype=int)

# Perform max pooling
for i in range(result_rows):
    for j in range(result_cols):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        result[i, j] = np.max(window)

# Print the result
print(result)