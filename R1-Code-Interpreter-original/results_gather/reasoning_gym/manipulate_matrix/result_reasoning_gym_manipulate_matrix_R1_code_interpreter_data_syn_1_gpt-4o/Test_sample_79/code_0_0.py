import numpy as np

# Initial matrix
matrix = np.array([
    [7, 3, 6, 6, 9, 2, 3],
    [6, 8, 3, 2, 7, 6, 6],
    [4, 6, 3, 5, 7, 0, 5],
    [3, 9, 7, 6, 1, 6, 8],
    [1, 4, 1, 5, 7, 3, 5]
])

# Step 1: Identity transformation (no change)
# Step 2: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, [1, 3], axis=0)

# Step 3: Crop the matrix to rows 2-3 and columns 6-7 (1-indexed)
matrix = matrix[1:3, 5:7]

# Step 4: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 5: Mirror the matrix along the counterdiagonal again
matrix = np.fliplr(np.flipud(matrix.T))

# Print the final matrix
print(matrix)