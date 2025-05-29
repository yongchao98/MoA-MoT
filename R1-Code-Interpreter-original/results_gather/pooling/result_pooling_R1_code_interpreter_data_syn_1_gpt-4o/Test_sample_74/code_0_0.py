import numpy as np

# Define the original 5x5 matrix
matrix = np.array([
    [-10, -7, -7, 4, 12],
    [-4, 14, 20, -5, 20],
    [3, -9, 3, -9, 5],
    [12, -1, -4, -2, 5],
    [6, -7, -6, 9, 10]
])

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the result matrix
result = np.zeros((output_size, output_size))

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the current window
        window_average = np.mean(window)
        # Store the result
        result[i, j] = window_average

# Print the resulting matrix
print(result)