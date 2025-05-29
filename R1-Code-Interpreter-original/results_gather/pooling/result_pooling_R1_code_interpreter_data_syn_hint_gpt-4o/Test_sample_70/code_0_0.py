# Define the original matrix
matrix = [
    [10, -5, 20, 11, 6],
    [5, 18, 7, 1, 16],
    [17, -2, 1, -2, 11],
    [7, 17, 14, 17, 9],
    [17, 15, 4, -1, -8]
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
        # Calculate the sum of the current window
        window_sum = 0
        for wi in range(window_size):
            for wj in range(window_size):
                window_sum += matrix[i + wi][j + wj]
        # Calculate the average and store it in the output matrix
        output_matrix[i][j] = window_sum / (window_size * window_size)

# Print the output matrix
print(output_matrix)