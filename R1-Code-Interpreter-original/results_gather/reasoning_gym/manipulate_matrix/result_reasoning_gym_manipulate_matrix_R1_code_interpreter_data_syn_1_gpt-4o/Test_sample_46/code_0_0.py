import numpy as np

# Initial matrix
matrix = np.array([
    [4, 4, 6],
    [0, 2, 6],
    [3, 7, 6],
    [8, 0, 3],
    [2, 5, 6]
])

# Step 1: Identity transformation (no change)
# Step 2: Remove every 3rd column (1-indexed)
matrix = np.delete(matrix, 2, axis=1)

# Step 3: Crop the matrix to rows 2-3 and columns 2-2 (1-indexed)
matrix = matrix[1:3, 1:2]

# Step 4: Rotate the matrix 90 degrees clockwise
matrix = np.rot90(matrix, -1)

# Step 5: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 6: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 7: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Print the final matrix
print(matrix)