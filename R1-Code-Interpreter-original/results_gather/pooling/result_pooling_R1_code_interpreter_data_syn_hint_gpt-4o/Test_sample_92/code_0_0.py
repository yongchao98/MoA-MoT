# Define the original matrix
matrix = [
    [13, -9, -3, 17, 2, -7],
    [16, 1, 18, 12, -8, 2],
    [17, -8, -1, 10, 2, 18],
    [-5, 4, -2, 12, -4, -9],
    [-2, -4, 6, -1, -5, 4],
    [3, 7, 11, -6, 0, -3]
]

# Dimensions of the sliding window
n = 3

# Calculate the dimensions of the output matrix
output_rows = len(matrix) - n + 1
output_cols = len(matrix[0]) - n + 1

# Initialize the output matrix
output_matrix = [[0] * output_cols for _ in range(output_rows)]

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Calculate the sum of the current 3x3 window
        window_sum = 0
        for k in range(n):
            for l in range(n):
                window_sum += matrix[i + k][j + l]
        # Calculate the average
        window_average = window_sum / (n * n)
        # Store the average in the output matrix
        output_matrix[i][j] = window_average

# Print the output matrix
print(output_matrix)