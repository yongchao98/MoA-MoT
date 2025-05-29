import numpy as np

# Original matrix
matrix = np.array([
    [-5, -10, 1],
    [1, 14, -10],
    [15, 19, 20]
])

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] - window_size + 1
output_cols = matrix.shape[1] - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        output_matrix[i, j] = np.max(window)

print(output_matrix)