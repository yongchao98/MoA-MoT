import numpy as np

# Define the original matrix
matrix = np.array([
    [0, 5, 3, 10, -7],
    [-10, -3, 11, -3, 19],
    [-3, 2, 5, -5, -1],
    [2, -9, 12, 8, 5],
    [3, 18, -6, 6, 15]
])

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_size, output_size))

# Perform the average pooling operation
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the current window
        window_average = np.mean(window)
        # Store the result in the output matrix
        output_matrix[i, j] = window_average

# Print the output matrix
print(output_matrix)