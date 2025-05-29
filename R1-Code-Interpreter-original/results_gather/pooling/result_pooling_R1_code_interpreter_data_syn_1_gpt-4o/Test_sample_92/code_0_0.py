# Define the input matrix
matrix = [
    [13, -9, -3, 17, 2, -7],
    [16, 1, 18, 12, -8, 2],
    [17, -8, -1, 10, 2, 18],
    [-5, 4, -2, 12, -4, -9],
    [-2, -4, 6, -1, -5, 4],
    [3, 7, 11, -6, 0, -3]
]

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_rows = len(matrix) - window_size + 1
output_cols = len(matrix[0]) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_cols for _ in range(output_rows)]

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Calculate the sum of the current 3x3 window
        window_sum = 0
        for wi in range(window_size):
            for wj in range(window_size):
                window_sum += matrix[i + wi][j + wj]
        # Calculate the average and store it in the output matrix
        output_matrix[i][j] = window_sum / (window_size * window_size)

# Print the output matrix
print(output_matrix)