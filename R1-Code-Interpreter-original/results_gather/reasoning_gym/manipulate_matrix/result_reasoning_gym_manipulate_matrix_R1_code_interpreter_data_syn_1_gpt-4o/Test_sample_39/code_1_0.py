import numpy as np

# Initial matrix
matrix = np.array([
    [9, 0, 4, 2, 3, 4, 1, 3, 8, 0],
    [9, 3, 8, 0, 8, 5, 7, 2, 5, 2],
    [0, 5, 9, 0, 6, 9, 3, 2, 5, 6],
    [7, 1, 0, 3, 9, 2, 2, 3, 8, 5],
    [4, 9, 0, 4, 5, 3, 6, 0, 1, 1],
    [3, 1, 4, 2, 6, 7, 7, 9, 9, 5],
    [4, 1, 3, 4, 9, 3, 2, 9, 7, 7]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 3: Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=3)

# Step 4: Crop the matrix to rows 7-7 and columns 1-4 (1-indexed)
matrix = matrix[6:7, 0:4]

# Step 5: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 6: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 7: Remove every 4-th row (1-indexed)
matrix = np.delete(matrix, np.arange(3, matrix.shape[0], 4), axis=0)

# Step 8: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[0], 2), axis=0)

# Step 9: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 10: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Convert the final matrix to the required format
result = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(f"<<<{result}>>>")