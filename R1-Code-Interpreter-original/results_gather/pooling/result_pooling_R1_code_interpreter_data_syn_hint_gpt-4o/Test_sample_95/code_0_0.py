import numpy as np

# Define the original matrix
matrix = [
    [7, 4, 11, 15, 18, 2],
    [14, 3, 1, 0, 19, 11],
    [5, 8, 14, -8, -10, 15],
    [-7, -9, 9, -8, 9, -6],
    [20, 5, 9, 1, 8, 19],
    [16, 12, -10, 2, 6, -1]
]

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_rows = len(matrix) - window_size + 1
output_cols = len(matrix[0]) - window_size + 1

# Initialize the output matrix
output_matrix = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current window
        window = [matrix[i+x][j+y] for x in range(window_size) for y in range(window_size)]
        # Calculate the average of the current window
        window_average = sum(window) / (window_size * window_size)
        # Store the average in the output matrix
        output_matrix[i][j] = window_average

# Print the output matrix
print(output_matrix)