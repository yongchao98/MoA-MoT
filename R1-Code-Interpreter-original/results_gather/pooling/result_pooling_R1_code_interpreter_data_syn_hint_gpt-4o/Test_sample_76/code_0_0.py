import numpy as np

# Define the original matrix
matrix = np.array([
    [15, 5, -4, 6, 7],
    [2, -4, -5, 4, -8],
    [-3, 18, -6, -9, -5],
    [-8, -6, -2, -6, 15],
    [4, 1, -2, -7, 13]
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
        # Extract the 3x3 sub-matrix
        sub_matrix = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average
        average = np.mean(sub_matrix)
        # Store the result
        new_matrix[i, j] = average

# Print the resulting matrix
print(new_matrix)