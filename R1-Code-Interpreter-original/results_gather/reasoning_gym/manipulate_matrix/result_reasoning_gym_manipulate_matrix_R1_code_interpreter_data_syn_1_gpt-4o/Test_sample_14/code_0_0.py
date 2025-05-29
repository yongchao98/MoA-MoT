import numpy as np

# Initial matrix
matrix = np.array([
    [0, 9, 9, 3, 8, 2, 9, 8],
    [0, 0, 1, 5, 0, 0, 2, 8],
    [4, 8, 5, 1, 6, 7, 5, 2]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 3: Rotate 90 degrees
matrix = np.rot90(matrix, -1)

# Step 4: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 5: Mirror along the diagonal
matrix = matrix.T

# Step 6: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 7: Crop to rows 3-3 and columns 3-3 (1-indexed)
matrix = matrix[2:3, 2:3]

# Step 8: Set elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Step 9: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Output the final matrix
print(matrix)