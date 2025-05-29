# Define the original matrix
matrix = [
    [2, -2, 10, 3, 3],
    [1, 3, 9, 15, 16],
    [20, 20, 6, 13, 15],
    [5, -3, 19, -1, 14],
    [7, 19, 4, 1, 16]
]

# Define the size of the sliding window
window_size = 3
N = len(matrix)

# Initialize the result matrix
result_matrix = []

# Perform the average pooling
for i in range(N - window_size + 1):
    row = []
    for j in range(N - window_size + 1):
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