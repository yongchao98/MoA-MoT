# Define the original 5x5 matrix
matrix = [
    [-8, -3, 6, -2, 17],
    [18, -8, 10, 8, -5],
    [10, 11, 19, 16, 3],
    [2, -9, 20, 5, 16],
    [14, -7, 13, 7, 15]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the resulting matrix after pooling
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