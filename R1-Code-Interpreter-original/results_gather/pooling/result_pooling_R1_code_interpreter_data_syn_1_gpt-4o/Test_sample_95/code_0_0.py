# Define the original 6x6 matrix
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

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform the average pooling operation
for i in range(result_size):
    for j in range(result_size):
        # Calculate the sum of the current 3x3 window
        window_sum = 0
        for wi in range(window_size):
            for wj in range(window_size):
                window_sum += matrix[i + wi][j + wj]
        # Calculate the average and store it in the result matrix
        result_matrix[i][j] = window_sum / (window_size * window_size)

# Print the result matrix
print(result_matrix)