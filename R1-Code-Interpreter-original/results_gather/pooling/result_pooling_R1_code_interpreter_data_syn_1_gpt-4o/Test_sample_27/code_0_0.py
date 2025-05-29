import numpy as np

# Original matrix
matrix = np.array([
    [11, 15, 18, 18],
    [20, -8, -4, 19],
    [14, -9, 11, -7],
    [-9, 2, -7, 2]
])

# Dimensions
N = 4
n = 2

# Resulting matrix dimensions
result_dim = N - n + 1
result_matrix = np.zeros((result_dim, result_dim), dtype=int)

# Perform max pooling
for i in range(result_dim):
    for j in range(result_dim):
        # Extract the 2x2 window
        window = matrix[i:i+n, j:j+n]
        # Find the maximum value in the window
        max_value = np.max(window)
        # Store the maximum value in the result matrix
        result_matrix[i, j] = max_value

# Print the result matrix
print(result_matrix)