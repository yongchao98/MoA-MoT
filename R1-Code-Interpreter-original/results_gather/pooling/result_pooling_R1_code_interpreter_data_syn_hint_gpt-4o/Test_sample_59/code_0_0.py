import numpy as np

# Define the original matrix
matrix = np.array([
    [11, 4, -7, -9],
    [3, -10, 2, 9],
    [-7, 5, -8, 13],
    [17, 1, 5, -9]
])

# Define the size of the sliding window
window_size = 2

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = np.zeros((result_size, result_size))

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Assign the average to the result matrix
        result_matrix[i, j] = window_average

# Print the result matrix
print(result_matrix)