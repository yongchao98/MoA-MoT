import numpy as np

# Original 4x4 matrix
matrix = np.array([
    [14, 10, 8, -3],
    [-4, 13, 8, 11],
    [-8, -4, -9, 1],
    [-9, 15, 8, 7]
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
        # Extract the 2x2 sub-matrix
        sub_matrix = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average
        average = np.mean(sub_matrix)
        # Store the average in the new matrix
        new_matrix[i, j] = average

# Print the resulting matrix
print(new_matrix)