# Define the original matrix
matrix = [
    [11, 4, -7, -9],
    [3, -10, 2, 9],
    [-7, 5, -8, 13],
    [17, 1, 5, -9]
]

# Dimensions of the original matrix
N = 4

# Size of the sliding window
n = 2

# Resulting matrix dimensions
result_dim = N - n + 1

# Initialize the resulting matrix
result_matrix = [[0] * result_dim for _ in range(result_dim)]

# Perform average pooling
for i in range(result_dim):
    for j in range(result_dim):
        # Calculate the sum of the current 2x2 window
        window_sum = (
            matrix[i][j] + matrix[i][j+1] +
            matrix[i+1][j] + matrix[i+1][j+1]
        )
        # Calculate the average
        result_matrix[i][j] = window_sum / (n * n)

# Print the resulting matrix
print(result_matrix)