import numpy as np

# Define the 6x6 matrix
matrix = np.array([
    [17, -10, 19, 2, -8, 11],
    [12, 1, 9, 19, 2, 11],
    [6, 3, 3, 12, 15, 15],
    [-9, 2, 14, 19, 6, -7],
    [14, 7, 2, 14, 5, 4],
    [19, -9, 13, 14, 14, 17]
])

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = matrix.shape[0] - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_size, output_size))

# Perform the average pooling operation
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