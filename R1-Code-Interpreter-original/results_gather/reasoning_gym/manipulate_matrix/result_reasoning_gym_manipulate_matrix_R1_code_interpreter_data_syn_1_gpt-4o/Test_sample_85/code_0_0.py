import numpy as np

# Initial matrix
matrix = np.array([
    [2, 3, 6, 2],
    [8, 6, 1, 8],
    [9, 5, 5, 7],
    [0, 2, 5, 6],
    [1, 2, 6, 6]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the counterdiagonal
matrix = np.fliplr(matrix).T

# Step 3: Horizontally mirror the matrix
matrix = np.flipud(matrix)

# Step 4: Mirror along the counterdiagonal again
matrix = np.fliplr(matrix).T

# Step 5: Horizontally mirror the matrix again
matrix = np.flipud(matrix)

# Step 6: Crop to rows 2-5 and columns 4-4 (1-indexed)
matrix = matrix[1:5, 3:4]

# Step 7: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, [1, 3], axis=0)

# Step 8: Mirror along the diagonal
matrix = matrix.T

# Print the final matrix
print(matrix)