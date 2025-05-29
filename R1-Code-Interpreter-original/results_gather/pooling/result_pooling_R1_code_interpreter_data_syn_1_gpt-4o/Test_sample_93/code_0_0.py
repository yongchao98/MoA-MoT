import numpy as np

# Define the original 6x6 matrix
matrix = np.array([
    [-7, 1, 0, -6, -9, -5],
    [8, 10, 8, 8, -4, -9],
    [8, 4, 8, 9, -1, -3],
    [2, 16, 5, -4, 8, -9],
    [-1, 17, 1, -4, -7, -2],
    [19, 18, 7, -9, 13, 10]
])

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
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