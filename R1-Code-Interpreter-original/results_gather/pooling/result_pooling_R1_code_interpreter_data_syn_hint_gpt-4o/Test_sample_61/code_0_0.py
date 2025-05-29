# Define the original matrix
matrix = [
    [-1, 11, -10, 8],
    [8, 13, -9, 13],
    [-8, 12, 5, 0],
    [-6, -6, 20, -9]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the resultant matrix
result_size = len(matrix) - window_size + 1

# Initialize the resultant matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average and store it in the resultant matrix
        result_matrix[i][j] = window_sum / (window_size * window_size)

# Print the resultant matrix
print(result_matrix)