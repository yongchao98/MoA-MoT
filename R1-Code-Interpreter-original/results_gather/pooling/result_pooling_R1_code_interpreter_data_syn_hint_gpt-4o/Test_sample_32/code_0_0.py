import numpy as np

# Original matrix
matrix = np.array([
    [-4, 11, 19, 19],
    [-8, -8, 1, -5],
    [-10, 4, -6, 13],
    [-4, 17, -4, -10]
])

# Dimensions of the sliding window
window_size = 2

# Dimensions of the original matrix
n = matrix.shape[0]

# Resulting matrix dimensions
result_size = n - window_size + 1

# Initialize the result matrix
result_matrix = np.zeros((result_size, result_size), dtype=int)

# Perform max pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current 2x2 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        max_value = np.max(window)
        # Store the maximum value in the result matrix
        result_matrix[i, j] = max_value

# Print the resulting matrix
print(result_matrix)