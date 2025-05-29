# Define the original matrix
matrix = [
    [-4, 10, 7, 14, 6, 5],
    [-5, 13, -3, -5, 1, 11],
    [-1, 5, 19, 5, 16, 15],
    [9, 0, -8, 19, 6, 6],
    [1, 12, 8, 3, -10, 15],
    [-5, 14, 18, 16, -6, 14]
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