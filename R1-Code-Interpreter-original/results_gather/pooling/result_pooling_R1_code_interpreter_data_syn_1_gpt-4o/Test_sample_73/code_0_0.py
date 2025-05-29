# Define the original 5x5 matrix
matrix = [
    [10, 8, -7, -3, 19],
    [-8, -5, 14, 11, -2],
    [14, -5, 9, 8, -9],
    [2, 10, -5, 7, 12],
    [18, 0, -9, -10, 1]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the result matrix
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