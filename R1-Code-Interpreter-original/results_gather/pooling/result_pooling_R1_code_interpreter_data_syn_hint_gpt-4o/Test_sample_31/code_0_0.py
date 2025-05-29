import numpy as np

# Original matrix
matrix = np.array([
    [5, 3, 2, 7],
    [-9, 13, 10, -3],
    [12, 3, -2, 16],
    [14, 4, 9, 18]
])

# Dimensions of the original matrix
N = matrix.shape[0]

# Size of the sliding window
n = 2

# Resulting matrix dimensions
result_dim = N - n + 1

# Initialize the resulting matrix
result_matrix = np.zeros((result_dim, result_dim), dtype=int)

# Perform max pooling
for i in range(result_dim):
    for j in range(result_dim):
        # Extract the 2x2 window
        window = matrix[i:i+n, j:j+n]
        # Find the maximum value in the window
        max_value = np.max(window)
        # Assign the maximum value to the resulting matrix
        result_matrix[i, j] = max_value

# Print the resulting matrix
print(result_matrix)