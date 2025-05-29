# Define the original matrix
matrix = [
    [0, 19, 1, 1, 17, 11],
    [13, 12, -4, 12, 17, -5],
    [18, -2, -8, -9, 20, -5],
    [-6, 2, -7, 0, 8, -10],
    [6, 2, 5, 0, 7, -8],
    [-1, 10, 8, 6, 0, 4]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the result matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform average pooling
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