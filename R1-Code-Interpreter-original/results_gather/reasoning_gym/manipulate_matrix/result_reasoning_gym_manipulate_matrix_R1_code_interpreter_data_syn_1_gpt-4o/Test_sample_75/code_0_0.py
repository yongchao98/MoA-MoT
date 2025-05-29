import numpy as np

# Initial matrix
matrix = np.array([
    [7, 4, 2, 6, 8, 1, 4, 4],
    [9, 0, 4, 5, 3, 2, 9, 3],
    [2, 8, 2, 3, 5, 1, 1, 2],
    [1, 7, 3, 2, 7, 4, 1, 3],
    [4, 6, 6, 3, 1, 9, 4, 6]
])

# Step 2: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, [1, 3], axis=0)

# Step 3: Remove every 3-th row (1-indexed)
matrix = np.delete(matrix, [2], axis=0)

# Step 4: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 5: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 6: Rotate the matrix 360 degrees (no change)
# No operation needed as 360-degree rotation results in the same matrix

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))