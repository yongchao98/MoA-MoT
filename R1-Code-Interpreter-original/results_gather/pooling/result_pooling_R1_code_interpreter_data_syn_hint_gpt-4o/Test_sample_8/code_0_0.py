import numpy as np

# Original matrix
matrix = np.array([
    [-9, 1, 18],
    [0, 19, -1],
    [10, 14, -1]
])

# Dimensions of the sliding window
window_size = 2

# Resulting matrix after max pooling
pooled_matrix = np.zeros((2, 2))

# Perform max pooling
for i in range(2):
    for j in range(2):
        # Extract the current 2x2 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        pooled_matrix[i, j] = np.max(window)

# Print the resulting pooled matrix
print(pooled_matrix)