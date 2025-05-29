import numpy as np

# Define the original matrix
matrix = np.array([
    [17, -10, 13],
    [7, -5, 9],
    [-9, 11, 17]
])

# Define the size of the sliding window
window_size = 2

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_size, output_size), dtype=int)

# Perform max pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Find the maximum value in the current window
        max_value = np.max(window)
        # Assign the maximum value to the output matrix
        output_matrix[i, j] = max_value

# Print the output matrix
print(output_matrix)