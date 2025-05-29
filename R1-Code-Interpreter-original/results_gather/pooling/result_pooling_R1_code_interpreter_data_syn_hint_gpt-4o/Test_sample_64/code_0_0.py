import numpy as np

# Original matrix
matrix = np.array([
    [7, 18, 18, -4],
    [12, 6, 14, -6],
    [1, -2, -8, 3],
    [20, 16, 9, -6]
])

# Dimensions
N = 4
n = 2

# Resultant matrix dimensions
result_dim = N - n + 1
result_matrix = np.zeros((result_dim, result_dim))

# Perform average pooling
for i in range(result_dim):
    for j in range(result_dim):
        window = matrix[i:i+n, j:j+n]
        average = np.mean(window)
        result_matrix[i, j] = average

print(result_matrix)