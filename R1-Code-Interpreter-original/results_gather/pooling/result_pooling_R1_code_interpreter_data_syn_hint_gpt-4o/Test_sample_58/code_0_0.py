# Define the original matrix
matrix = [
    [-10, 9, -6, -2],
    [-9, 9, 7, 15],
    [12, -3, -10, 7],
    [5, 3, 17, 9]
]

# Dimensions of the original matrix
N = 4

# Dimensions of the sliding window
n = 2

# Resultant matrix dimensions
result_dim = N - n + 1

# Initialize the result matrix
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
        window_avg = window_sum / (n * n)
        # Store the average in the result matrix
        result_matrix[i][j] = window_avg

# Print the result matrix
print(result_matrix)