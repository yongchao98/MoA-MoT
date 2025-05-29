# Original matrix
matrix = [
    [1, 4, 11, 3, 12],
    [-9, 13, 11, 1, 9],
    [16, 12, -9, 7, 4],
    [8, 0, 14, -5, 6],
    [-4, 0, -7, -1, 1]
]

# Dimensions
N = 5
n = 3

# Resulting matrix dimensions
result_dim = N - n + 1

# Initialize the resulting matrix
result_matrix = [[0] * result_dim for _ in range(result_dim)]

# Perform max pooling
for i in range(result_dim):
    for j in range(result_dim):
        # Extract the 3x3 window
        window = [matrix[x][j:j+n] for x in range(i, i+n)]
        # Find the maximum value in the window
        max_value = max(max(row) for row in window)
        # Assign the max value to the resulting matrix
        result_matrix[i][j] = max_value

# Print the resulting matrix
print(result_matrix)