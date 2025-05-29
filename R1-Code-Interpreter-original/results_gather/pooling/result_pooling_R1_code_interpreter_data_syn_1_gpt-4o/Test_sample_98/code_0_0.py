# Define the original 6x6 matrix
matrix = [
    [-10, 19, -7, 11, -3, -3],
    [20, -9, -5, 0, 15, -1],
    [12, 4, 16, -8, 16, -1],
    [2, 7, 20, -8, -10, 6],
    [19, 4, -8, 5, 15, 3],
    [6, 19, 2, -5, -6, 20]
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