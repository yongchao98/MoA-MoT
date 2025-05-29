import numpy as np

# Initial matrix
matrix = np.array([
    [4, 8, 3, 7, 7, 8],
    [7, 8, 5, 6, 3, 1],
    [8, 7, 4, 3, 7, 9],
    [0, 4, 1, 2, 8, 2],
    [1, 2, 5, 5, 9, 5],
    [9, 3, 3, 6, 7, 7],
    [0, 2, 0, 5, 4, 5],
    [4, 0, 2, 6, 1, 4]
])

# Identity transformation (no change)
# Rotate the matrix 180 degrees
matrix = np.rot90(matrix, 2)

# Map each occurrence of 8 to 6
matrix[matrix == 8] = 6

# Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Map each occurrence of 6 to 5
matrix[matrix == 6] = 5

# Vertically mirror the matrix
matrix = np.flipud(matrix)

# Vertically mirror the matrix again
matrix = np.flipud(matrix)

# Mirror the matrix along the diagonal
matrix = matrix.T

# Mirror the matrix along the diagonal again
matrix = matrix.T

# Set all elements divisible by 1 to zero
matrix = np.zeros_like(matrix)

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))