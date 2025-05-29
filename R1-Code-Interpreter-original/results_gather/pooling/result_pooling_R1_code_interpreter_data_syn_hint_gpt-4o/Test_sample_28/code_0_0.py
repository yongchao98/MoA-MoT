# Original matrix
matrix = [
    [9, -9, 4, -3],
    [-1, -8, 15, 7],
    [11, 13, -9, 17],
    [20, 15, 20, -7]
]

# Dimensions of the sliding window
window_size = 2

# Dimensions of the original matrix
N = len(matrix)

# Resultant matrix dimensions
result_size = N - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform max pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the window
        result_matrix[i][j] = max(window)

# Print the result matrix
print(result_matrix)