# Define the original matrix
matrix = [
    [3, 13, 4, 14],
    [17, -9, -5, 10],
    [8, 0, 7, 13],
    [13, -6, 1, 8]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average and store it in the result matrix
        result_matrix[i][j] = window_sum / (window_size * window_size)

# Print the result matrix
print(result_matrix)