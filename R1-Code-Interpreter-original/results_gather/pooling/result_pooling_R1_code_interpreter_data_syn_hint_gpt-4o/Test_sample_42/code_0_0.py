import numpy as np

# Original matrix
matrix = np.array([
    [2, -7, -6, 6, 0],
    [-8, 13, -2, 8, -7],
    [20, 9, -8, 17, 3],
    [11, -8, 19, 14, 3],
    [-3, 15, 0, 8, 9]
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
        # Extract the 3x3 window
        window = matrix[i:i+n, j:j+n]
        # Find the maximum value in the window
        max_value = np.max(window)
        # Assign the max value to the result matrix
        result_matrix[i, j] = max_value

# Print the resulting matrix
print(result_matrix)