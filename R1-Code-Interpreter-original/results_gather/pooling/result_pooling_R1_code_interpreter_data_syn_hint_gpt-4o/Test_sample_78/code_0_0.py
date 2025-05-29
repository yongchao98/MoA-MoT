import numpy as np

# Define the 5x5 matrix
matrix = np.array([
    [-2, 16, 3, -4, 13],
    [14, -9, -6, -5, 1],
    [18, -10, -4, 15, 20],
    [7, 4, 18, 20, 14],
    [19, 14, -5, -1, 16]
])

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
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