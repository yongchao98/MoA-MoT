import numpy as np

# Original matrix
matrix = np.array([
    [8, 10, -6, -7],
    [0, 10, 2, 19],
    [4, 3, -5, -10],
    [2, 2, -1, -9]
])

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the new matrix
new_matrix_size = len(matrix) - window_size + 1

# Initialize the new matrix
new_matrix = np.zeros((new_matrix_size, new_matrix_size))

# Perform average pooling
for i in range(new_matrix_size):
    for j in range(new_matrix_size):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Assign the average to the new matrix
        new_matrix[i, j] = window_average

# Print the resulting matrix
print(new_matrix)