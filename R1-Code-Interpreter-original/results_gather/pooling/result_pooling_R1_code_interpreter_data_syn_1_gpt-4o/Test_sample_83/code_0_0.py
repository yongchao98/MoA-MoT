# Define the original 5x5 matrix
matrix = [
    [-5, -2, -1, 9, 11],
    [9, 6, 12, 16, 16],
    [11, 14, 18, 0, 16],
    [20, 14, 15, -7, 3],
    [5, 6, -1, 1, 7]
]

# Dimensions of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform the average pooling operation
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