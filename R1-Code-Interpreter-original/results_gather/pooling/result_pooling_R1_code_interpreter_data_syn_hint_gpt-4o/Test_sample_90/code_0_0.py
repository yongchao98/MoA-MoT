import numpy as np

# Define the 6x6 matrix
matrix = np.array([
    [0, 19, 1, 1, 17, 11],
    [13, 12, -4, 12, 17, -5],
    [18, -2, -8, -9, 20, -5],
    [-6, 2, -7, 0, 8, -10],
    [6, 2, 5, 0, 7, -8],
    [-1, 10, 8, 6, 0, 4]
])

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_size, output_size))

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the result in the output matrix
        output_matrix[i, j] = window_average

# Print the output matrix
print(output_matrix)