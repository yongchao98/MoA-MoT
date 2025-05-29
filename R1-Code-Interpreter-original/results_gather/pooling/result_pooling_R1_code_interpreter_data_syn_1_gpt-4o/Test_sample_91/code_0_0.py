import numpy as np

# Define the 6x6 matrix
matrix = np.array([
    [5, 8, 3, 6, 3, 16],
    [2, 8, 7, -2, 10, 19],
    [-1, 8, -6, 8, 15, -3],
    [0, -3, 18, 2, 12, -4],
    [19, 10, -6, -10, 1, -5],
    [16, 8, 20, 1, -10, 17]
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