import numpy as np

# Initial matrix
matrix = np.array([
    [9, 6, 2, 9, 8],
    [2, 5, 4, 2, 6],
    [0, 9, 6, 2, 0]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 3: Remove every 4-th row (1-indexed)
# No rows to remove since there are only 3 rows

# Step 4: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 5: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 6: Mirror along the diagonal
matrix = matrix.T

# Step 7: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 8: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 9: Remove every 4-th column (1-indexed)
matrix = np.delete(matrix, 3, axis=1)

# Print the final matrix
print(matrix)