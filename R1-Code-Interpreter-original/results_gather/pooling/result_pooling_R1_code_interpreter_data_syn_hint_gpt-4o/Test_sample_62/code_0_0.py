# Define the original matrix
matrix = [
    [9, 11, 15, -9],
    [8, -3, 1, 5],
    [4, -9, 7, 11],
    [7, -8, 0, 6]
]

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the resulting matrix
result_rows = len(matrix) - window_size + 1
result_cols = len(matrix[0]) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_cols for _ in range(result_rows)]

# Perform average pooling
for i in range(result_rows):
    for j in range(result_cols):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average
        result_matrix[i][j] = window_sum / (window_size * window_size)

# Print the result matrix
print(result_matrix)