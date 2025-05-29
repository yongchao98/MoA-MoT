import numpy as np

# Original matrix
matrix = np.array([
    [17, 10, 6, 10],
    [0, -5, 16, 16],
    [0, 8, 0, 19],
    [-2, 17, -8, 11]
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
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        output_matrix[i, j] = np.max(window)

# Print the resulting matrix
print(output_matrix)