import numpy as np

# Original matrix
matrix = np.array([
    [0, 9, -1, -10, 20],
    [4, 2, -7, -3, 6],
    [-5, -10, 5, 14, 3],
    [7, -4, 7, -1, 6],
    [7, 5, 2, -9, 17]
])

# Dimensions
N = 5
n = 3

# Resulting matrix dimensions
result_dim = N - n + 1
result_matrix = np.zeros((result_dim, result_dim))

# Perform max pooling
for i in range(result_dim):
    for j in range(result_dim):
        # Extract the current 3x3 window
        window = matrix[i:i+n, j:j+n]
        # Find the maximum value in the current window
        max_value = np.max(window)
        # Assign the maximum value to the result matrix
        result_matrix[i, j] = max_value

# Print the resulting matrix
print(result_matrix)