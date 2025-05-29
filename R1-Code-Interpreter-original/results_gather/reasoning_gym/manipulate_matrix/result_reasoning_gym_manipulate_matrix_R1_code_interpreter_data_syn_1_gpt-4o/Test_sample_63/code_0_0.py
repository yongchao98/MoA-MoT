import numpy as np

# Initial matrix
matrix = np.array([
    [8, 2, 7, 0, 2, 0, 9, 1, 4, 3],
    [0, 9, 8, 9, 3, 1, 5, 4, 6, 5],
    [2, 2, 2, 6, 5, 9, 2, 9, 6, 9],
    [4, 8, 2, 1, 1, 5, 7, 4, 0, 0],
    [5, 8, 8, 9, 5, 4, 5, 3, 6, 9],
    [6, 1, 2, 1, 9, 0, 0, 3, 6, 3],
    [8, 7, 9, 2, 6, 7, 2, 9, 9, 1],
    [7, 7, 3, 9, 5, 4, 6, 8, 3, 4],
    [2, 3, 1, 2, 8, 9, 3, 1, 0, 3],
    [7, 0, 5, 4, 2, 6, 3, 3, 8, 0]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 3: Crop to rows 3-3 and columns 5-8 (1-indexed)
matrix = matrix[4:5, 2:6]

# Step 4: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 5: Mirror along the counterdiagonal (reverse rows and columns)
matrix = np.fliplr(np.flipud(matrix))

# Step 6: Rotate 360 degrees (no change)
# Step 7: Crop to rows 1-1 and columns 3-3 (1-indexed)
matrix = matrix[0:1, 2:3]

# Step 8: Rotate 180 degrees (reverse rows and columns)
matrix = np.fliplr(np.flipud(matrix))

# Step 9: Mirror along the counterdiagonal (reverse rows and columns)
matrix = np.fliplr(np.flipud(matrix))

# Step 10: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 11: Rotate 180 degrees (reverse rows and columns)
matrix = np.fliplr(np.flipud(matrix))

# Output the final matrix
print(matrix)