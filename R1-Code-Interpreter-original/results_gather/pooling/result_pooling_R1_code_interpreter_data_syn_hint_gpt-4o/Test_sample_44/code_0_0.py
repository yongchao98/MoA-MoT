import numpy as np

# Define the 5x5 matrix
matrix = np.array([
    [20, 9, 7, 20, 7],
    [3, 19, -4, 9, 16],
    [7, 9, 7, -1, 10],
    [-6, 16, -1, -1, -3],
    [1, 7, -1, 12, 18]
])

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the result matrix
result = np.zeros((output_size, output_size), dtype=int)

# Perform max pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        result[i, j] = np.max(window)

# Print the result matrix
print(result)