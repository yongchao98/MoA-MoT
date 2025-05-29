import numpy as np

# Initial matrix
matrix = np.array([
    [0, 3, 4, 3, 2, 5, 4, 7, 0, 8],
    [2, 0, 3, 9, 7, 0, 2, 9, 5, 1],
    [4, 5, 2, 8, 8, 1, 6, 0, 5, 5],
    [3, 2, 4, 2, 8, 4, 4, 7, 2, 9],
    [3, 7, 6, 8, 1, 2, 9, 3, 9, 9],
    [0, 1, 7, 9, 4, 1, 3, 3, 7, 7],
    [1, 7, 5, 9, 7, 4, 8, 9, 3, 1],
    [0, 2, 9, 6, 2, 7, 1, 6, 3, 9],
    [0, 6, 7, 8, 3, 3, 7, 3, 2, 2]
])

# Remove every 3rd column (1-indexed)
matrix = np.delete(matrix, [2, 5, 8], axis=1)

# Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Crop the matrix to rows 1-2 and columns 3-3 (1-indexed)
matrix = matrix[0:2, 2:3]

# Remove every 2nd row (1-indexed)
matrix = np.delete(matrix, [1], axis=0)

# Print the final matrix
print(matrix)