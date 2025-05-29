import numpy as np

# Define the matrix
matrix = np.array([
    [8, 7, 9],
    [2, 7, 4],
    [9, 1, 7],
    [6, 8, 8],
    [3, 7, 1],
    [0, 4, 6]
])

# Identity transformation (no change)
identity_matrix = matrix.copy()

# Mirror along the counterdiagonal
# For a matrix of size m x n, the counterdiagonal mirroring involves swapping (i, j) with (n-j-1, m-i-1)
# However, since the matrix is not square, we need to handle it carefully
# We will create a new matrix with the same dimensions and fill it accordingly

rows, cols = matrix.shape
mirrored_matrix = np.zeros((cols, rows), dtype=int)

for i in range(rows):
    for j in range(cols):
        mirrored_matrix[j, rows - i - 1] = matrix[i, j]

# Print the final mirrored matrix
for row in mirrored_matrix:
    print(' '.join(map(str, row)))