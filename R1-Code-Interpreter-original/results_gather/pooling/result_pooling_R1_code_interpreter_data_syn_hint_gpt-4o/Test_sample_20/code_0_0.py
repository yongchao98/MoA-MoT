import numpy as np

# Original matrix
matrix = np.array([
    [1, -3, -9, 4],
    [-3, 15, -5, 2],
    [-5, 5, 2, 10],
    [-4, 6, 12, 13]
])

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] - window_size + 1
output_cols = matrix.shape[1] - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols), dtype=int)

# Perform max pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 2x2 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the window
        output_matrix[i, j] = np.max(window)

# Print the resulting matrix
print(output_matrix)