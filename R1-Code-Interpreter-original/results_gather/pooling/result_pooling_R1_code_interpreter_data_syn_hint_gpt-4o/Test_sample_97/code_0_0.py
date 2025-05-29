import numpy as np

# Define the 6x6 matrix
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

# Calculate the dimensions of the output matrix
output_rows = matrix.shape[0] - window_size + 1
output_cols = matrix.shape[1] - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current window
        window = matrix[i:i+window_size, j:j+window_size]
        # Calculate the average of the window
        window_average = np.mean(window)
        # Store the average in the output matrix
        output_matrix[i, j] = window_average

# Print the resulting matrix
print(output_matrix)