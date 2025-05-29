import numpy as np

# Define the original 6x6 matrix
matrix = np.array([
    [-8, 1, 1, -8, 5, -7],
    [4, 10, 19, 11, 17, 15],
    [5, 10, 5, -10, 13, -2],
    [-4, -10, 2, -9, 15, -8],
    [12, 11, 11, -7, 11, 14],
    [13, -9, 10, 5, -6, 10]
])

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = matrix.shape[0] - window_size + 1

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