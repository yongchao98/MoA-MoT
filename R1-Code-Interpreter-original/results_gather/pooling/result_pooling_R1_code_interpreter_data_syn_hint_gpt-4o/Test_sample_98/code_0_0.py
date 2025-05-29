import numpy as np

# Define the 6x6 matrix
matrix = np.array([
    [-10, 19, -7, 11, -3, -3],
    [20, -9, -5, 0, 15, -1],
    [12, 4, 16, -8, 16, -1],
    [2, 7, 20, -8, -10, 6],
    [19, 4, -8, 5, 15, 3],
    [6, 19, 2, -5, -6, 20]
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