# Define the original matrix
matrix = [
    [8, 8, 9, -3],
    [14, 10, 12, -9],
    [18, 11, 8, -10],
    [18, 18, 6, 19]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the new matrix
new_size = len(matrix) - window_size + 1

# Initialize the new matrix to store the averages
new_matrix = [[0] * new_size for _ in range(new_size)]

# Perform average pooling
for i in range(new_size):
    for j in range(new_size):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average and store it in the new matrix
        new_matrix[i][j] = window_sum / (window_size * window_size)

# Print the resulting matrix
print(new_matrix)