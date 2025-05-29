import numpy as np

# Initial matrix
matrix = np.array([
    [5, 3, 3],
    [7, 5, 9],
    [3, 6, 8],
    [8, 3, 9],
    [9, 9, 7]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 3: Rotate 90 degrees
matrix = np.rot90(matrix, -1)

# Step 4: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 5: Set elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Step 6: Crop to rows 2-2 and columns 2-3 (1-indexed)
matrix = matrix[1:2, 1:3]

# Step 7: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 8: Mirror along the diagonal
matrix = matrix.T

# Step 9: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 10: Mirror along the diagonal
matrix = matrix.T

# Step 11: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Print the final matrix
print(matrix)