import numpy as np

# Initial matrix
matrix = np.array([
    [8, 5, 2, 7, 0, 2, 5, 8],
    [4, 9, 8, 4, 1, 1, 0, 4],
    [5, 1, 7, 2, 6, 8, 5, 7],
    [5, 7, 5, 8, 5, 6, 4, 4],
    [3, 9, 0, 2, 9, 2, 1, 1]
])

# Step 1: Identity transformation (no change)

# Step 2: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 3: Mirror along the counterdiagonal again (revert to original)
matrix = np.fliplr(np.flipud(matrix.T))

# Step 4: Map each occurrence of 5 to 6
matrix[matrix == 5] = 6

# Step 5: Mirror along the diagonal
matrix = matrix.T

# Step 6: Remove every 5-th row (1-indexed)
matrix = np.delete(matrix, 4, axis=0)

# Step 7: Crop to rows 4-4 and columns 4-5 (1-indexed)
matrix = matrix[3:4, 3:5]

# Step 8: Map each occurrence of 5 to 7
matrix[matrix == 5] = 7

# Step 9: Crop to rows 1-1 and columns 1-2 (1-indexed)
matrix = matrix[0:1, 0:2]

# Step 10: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 11: Rotate the matrix 90 degrees
matrix = np.rot90(matrix, -1)

# Print the final matrix
print(matrix)