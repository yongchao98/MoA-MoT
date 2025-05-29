import numpy as np

# Initial matrix
matrix = np.array([
    [4, 7, 7, 0, 5, 3, 9, 5],
    [3, 5, 8, 1, 8, 3, 5, 0],
    [6, 0, 5, 2, 0, 7, 6, 5],
    [7, 6, 6, 1, 0, 8, 8, 8],
    [8, 2, 8, 2, 9, 4, 1, 2],
    [2, 1, 7, 6, 1, 4, 6, 6],
    [1, 5, 8, 2, 2, 5, 1, 4],
    [6, 8, 9, 6, 2, 6, 4, 6],
    [3, 3, 2, 2, 8, 5, 8, 7],
    [9, 3, 8, 2, 9, 4, 0, 5]
])

# Step 3: Mirror along the diagonal
matrix = matrix.T

# Step 4: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix))

# Step 5: Rotate 90 degrees
matrix = np.rot90(matrix, -1)

# Step 6: Vertically mirror
matrix = np.flipud(matrix)

# Step 8: Vertically mirror again
matrix = np.flipud(matrix)

# Step 9: Mirror along the diagonal again
matrix = matrix.T

# Step 10: Remove every 5th row (1-indexed)
# After transformations, the matrix should have 10 rows, so we remove the 5th and 10th rows
matrix = np.delete(matrix, [4, 9], axis=0)

# Print the final matrix
print(matrix)