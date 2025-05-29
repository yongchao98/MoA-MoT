import numpy as np

# Original matrix
matrix = np.array([
    [3, 13, 4, 14],
    [17, -9, -5, 10],
    [8, 0, 7, 13],
    [13, -6, 1, 8]
])

# Size of the sliding window
window_size = 2

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the resulting matrix
result_matrix = np.zeros((result_size, result_size))

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the average in the result matrix
        result_matrix[i, j] = window_average

# Print the resulting matrix
print(result_matrix)