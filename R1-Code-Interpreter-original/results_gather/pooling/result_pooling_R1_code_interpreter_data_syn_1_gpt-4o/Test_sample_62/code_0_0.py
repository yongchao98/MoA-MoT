# Define the original matrix
matrix = [
    [9, 11, 15, -9],
    [8, -3, 1, 5],
    [4, -9, 7, 11],
    [7, -8, 0, 6]
]

# Dimensions of the sliding window
window_size = 2

# Calculate the dimensions of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average
        window_average = window_sum / (window_size * window_size)
        # Store the average in the output matrix
        output_matrix[i][j] = window_average

# Print the output matrix
print(output_matrix)