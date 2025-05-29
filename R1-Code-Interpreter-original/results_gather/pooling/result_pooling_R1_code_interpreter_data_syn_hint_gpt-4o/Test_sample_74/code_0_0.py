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

# Calculate the dimensions of the new matrix
new_matrix_size = len(matrix) - window_size + 1

# Initialize the new matrix to store the averages
new_matrix = np.zeros((new_matrix_size, new_matrix_size))

# Perform the average pooling
for i in range(new_matrix_size):
    for j in range(new_matrix_size):
        # Extract the current 3x3 window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the average in the new matrix
        new_matrix[i, j] = window_average

# Print the resulting matrix
print(new_matrix)