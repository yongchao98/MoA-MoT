# Define the original 6x6 matrix
matrix = [
    [2, 4, 11, -2, -10, 8],
    [0, 5, 3, 15, -7, 10],
    [20, 2, 0, 4, 0, 6],
    [14, 9, -5, 10, 6, 17],
    [-7, 17, 4, -3, 20, 0],
    [8, 9, 13, -10, -2, 11]
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