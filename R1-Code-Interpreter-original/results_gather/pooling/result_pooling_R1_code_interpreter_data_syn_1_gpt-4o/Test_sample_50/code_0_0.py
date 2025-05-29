# Define the original matrix
matrix = [
    [2, -7, -10, 1],
    [14, 5, 18, -5],
    [7, 7, 18, -9],
    [9, 9, 10, -5]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the new matrix
new_size = len(matrix) - window_size + 1

# Initialize the new matrix
new_matrix = [[0] * new_size for _ in range(new_size)]

# Perform average pooling
for i in range(new_size):
    for j in range(new_size):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average
        new_matrix[i][j] = window_sum / (window_size * window_size)

# Print the new matrix
print(new_matrix)