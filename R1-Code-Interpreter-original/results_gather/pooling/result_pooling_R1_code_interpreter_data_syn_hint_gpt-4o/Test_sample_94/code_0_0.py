# Define the original matrix
matrix = [
    [17, -10, 19, 2, -8, 11],
    [12, 1, 9, 19, 2, 11],
    [6, 3, 3, 12, 15, 15],
    [-9, 2, 14, 19, 6, -7],
    [14, 7, 2, 14, 5, 4],
    [19, -9, 13, 14, 14, 17]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Calculate the sum of the current 3x3 window
        window_sum = 0
        for wi in range(window_size):
            for wj in range(window_size):
                window_sum += matrix[i + wi][j + wj]
        # Calculate the average and store it in the output matrix
        output_matrix[i][j] = window_sum / (window_size * window_size)

# Print the output matrix
print(output_matrix)