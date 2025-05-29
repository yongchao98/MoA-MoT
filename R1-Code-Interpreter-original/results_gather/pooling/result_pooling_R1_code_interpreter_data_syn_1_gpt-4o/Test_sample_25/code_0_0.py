# Define the original matrix
matrix = [
    [16, 4, 15, 5],
    [16, 15, 4, 6],
    [18, 1, -5, -2],
    [8, 20, -10, 14]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform max pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the current window
        result_matrix[i][j] = max(window)

# Print the result matrix
print(result_matrix)