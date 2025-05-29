import numpy as np

# Initial matrix
matrix = np.array([
    [9, 0, 3, 0, 6, 5, 1],
    [6, 4, 5, 1, 9, 8, 2],
    [1, 4, 5, 9, 0, 2, 8],
    [1, 0, 5, 4, 8, 0, 6],
    [2, 9, 4, 4, 6, 5, 3],
    [3, 7, 8, 2, 9, 2, 9],
    [6, 8, 6, 3, 3, 1, 9],
    [1, 7, 4, 6, 1, 9, 4],
    [3, 2, 6, 5, 7, 1, 8],
    [7, 6, 7, 5, 2, 7, 4]
])

# Identity transformation (no change)
matrix_identity = matrix.copy()

# Vertically mirror the matrix
matrix_vertical_mirror = np.fliplr(matrix_identity)

# Vertically mirror the matrix again (revert to original)
matrix_reverted = np.fliplr(matrix_vertical_mirror)

# Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix_rotated_270 = np.rot90(matrix_reverted, k=3)

# Print the final matrix
for row in matrix_rotated_270:
    print(' '.join(map(str, row)))