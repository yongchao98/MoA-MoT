import numpy as np

# Define the input matrix
matrix = np.array([
    [8, 10, -6, -7],
    [0, 10, 2, 19],
    [4, 3, -5, -10],
    [2, 2, -1, -9]
])

# Define the size of the sliding window
window_size = 2

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
        # Calculate the average of the current window
        window_average = np.mean(window)
        # Store the result in the output matrix
        output_matrix[i, j] = window_average

# Print the output matrix
print(output_matrix)