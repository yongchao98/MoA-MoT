import numpy as np

# Initial matrix
matrix = np.array([
    [5, 2, 0, 8, 0],
    [5, 7, 2, 3, 1],
    [2, 9, 2, 7, 6],
    [5, 6, 2, 4, 2],
    [1, 7, 5, 0, 2],
    [3, 7, 3, 3, 2],
    [0, 3, 8, 7, 8],
    [0, 7, 6, 8, 6]
])

# Step 2: Map each occurrence of 7 to 2
matrix[matrix == 7] = 2

# Step 3: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 4: Horizontally mirror the matrix (reverse each row)
matrix = np.fliplr(matrix)

# Step 5: Mirror the matrix along the counterdiagonal (reverse both rows and columns)
matrix = np.flipud(np.fliplr(matrix))

# Step 6: Mirror the matrix along the diagonal (transpose again)
matrix = matrix.T

# Step 7: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[0], 2), axis=0)

# Step 8: Remove every 2-nd row (1-indexed) again
matrix = np.delete(matrix, np.arange(1, matrix.shape[0], 2), axis=0)

# Step 9: Map each occurrence of 1 to 6
matrix[matrix == 1] = 6

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))