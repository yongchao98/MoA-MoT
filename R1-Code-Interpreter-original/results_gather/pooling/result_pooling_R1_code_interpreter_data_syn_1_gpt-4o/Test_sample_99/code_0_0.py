# Define the original matrix
matrix = [
    [3, -10, -3, -1, -10, 8],
    [9, 12, -10, 13, -2, 4],
    [16, 7, 20, 1, 6, 7],
    [6, 11, 11, 3, 19, 1],
    [-2, 1, 11, -4, 8, -4],
    [20, -1, -1, -10, -7, 12]
]

# Define the size of the sliding window
window_size = 3
n = len(matrix)  # Size of the original matrix

# Initialize the result matrix
result_matrix = []

# Perform average pooling
for i in range(n - window_size + 1):
    row = []
    for j in range(n - window_size + 1):
        # Calculate the sum of the current 3x3 window
        window_sum = 0
        for wi in range(window_size):
            for wj in range(window_size):
                window_sum += matrix[i + wi][j + wj]
        # Calculate the average
        window_average = window_sum / (window_size * window_size)
        row.append(window_average)
    result_matrix.append(row)

# Print the result matrix
print(result_matrix)