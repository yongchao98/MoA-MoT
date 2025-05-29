import numpy as np

# Initial matrix
matrix = np.array([
    [0, 1, 2, 7, 7, 3],
    [4, 6, 9, 6, 9, 1],
    [7, 1, 5, 8, 5, 7],
    [0, 8, 8, 1, 6, 7],
    [2, 8, 5, 2, 7, 8],
    [1, 8, 7, 9, 3, 9],
    [6, 8, 4, 5, 6, 0],
    [6, 0, 6, 3, 5, 9],
    [1, 0, 6, 5, 4, 3],
    [5, 7, 6, 2, 4, 0]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 3: Rotate 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=1)

# Step 4: Rotate 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=1)

# Step 5: Remove every 6-th column (1-indexed)
matrix = np.delete(matrix, 5, axis=1)

# Step 6: Rotate 90 degrees (90 degrees clockwise)
matrix = np.rot90(matrix, k=-1)

# Step 7: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 8: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Convert the final matrix to the required format
result = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(f"<<<{result}>>>")