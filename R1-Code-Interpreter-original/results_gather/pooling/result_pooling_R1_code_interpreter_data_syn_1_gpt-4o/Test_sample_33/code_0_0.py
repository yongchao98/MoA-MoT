# Define the original matrix
matrix = [
    [17, 10, 6, 10],
    [0, -5, 16, 16],
    [0, 8, 0, 19],
    [-2, 17, -8, 11]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the resulting matrix
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

# Print the resulting matrix
print(result_matrix)