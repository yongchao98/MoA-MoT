# Define the original 6x6 matrix
matrix = [
    [5, 8, 3, 6, 3, 16],
    [2, 8, 7, -2, 10, 19],
    [-1, 8, -6, 8, 15, -3],
    [0, -3, 18, 2, 12, -4],
    [19, 10, -6, -10, 1, -5],
    [16, 8, 20, 1, -10, 17]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the resulting matrix
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