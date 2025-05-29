# Define the original matrix
matrix = [
    [19, 14, 3, 6],
    [5, 12, 7, -2],
    [-8, 13, -2, 20],
    [11, 10, -4, -10]
]

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the new matrix
new_matrix_size = len(matrix) - window_size + 1

# Initialize the new matrix
new_matrix = [[0] * new_matrix_size for _ in range(new_matrix_size)]

# Perform average pooling
for i in range(new_matrix_size):
    for j in range(new_matrix_size):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average
        new_matrix[i][j] = window_sum / (window_size * window_size)

# Print the resulting matrix
print(new_matrix)