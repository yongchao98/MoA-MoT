import numpy as np

# Initial matrix
matrix = np.array([
    [3, 3, 3, 1, 4, 7, 8],
    [9, 3, 2, 8, 8, 0, 4],
    [0, 7, 2, 3, 8, 6, 2],
    [5, 5, 1, 1, 8, 4, 8],
    [7, 3, 3, 3, 1, 1, 3],
    [8, 7, 9, 8, 1, 1, 3]
])

# Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, [1, 3, 5], axis=0)

# Map each occurrence of 1 to 8
matrix[matrix == 1] = 8

# Remove every 4-th column (1-indexed)
matrix = np.delete(matrix, 3, axis=1)

# Mirror the matrix along the counterdiagonal
# This involves reversing the order of rows and columns, then transposing
matrix = np.transpose(np.fliplr(np.flipud(matrix)))

# The matrix is already reduced to a single element
# Final matrix
print(matrix)