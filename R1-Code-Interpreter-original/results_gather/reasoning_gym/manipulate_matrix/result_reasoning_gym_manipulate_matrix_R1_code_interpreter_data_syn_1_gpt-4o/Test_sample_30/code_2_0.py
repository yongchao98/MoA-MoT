import numpy as np

# Define the original matrix
matrix = np.array([
    [3, 7, 4],
    [5, 2, 8],
    [2, 2, 0],
    [4, 2, 1]
])

# Identity transformation (no change)
identity_matrix = matrix.copy()

# Mirror along the counterdiagonal
# For a matrix of size m x n, the counterdiagonal mirroring involves swapping
# element (i, j) with element (n-j-1, m-i-1)

# Create a new matrix to store the mirrored result
mirrored_matrix = np.zeros_like(identity_matrix)

# Get the dimensions of the matrix
m, n = identity_matrix.shape

# Perform the mirroring
for i in range(m):
    for j in range(n):
        mirrored_matrix[i, j] = identity_matrix[m-j-1, n-i-1]

# Print the final mirrored matrix
print(mirrored_matrix)