import numpy as np

# Initial matrix
matrix = np.array([
    [6, 7, 4, 8, 1],
    [8, 6, 8, 5, 0],
    [1, 8, 1, 2, 7],
    [0, 8, 6, 9, 3],
    [4, 4, 5, 5, 8],
    [1, 6, 8, 7, 4],
    [2, 3, 4, 6, 7],
    [3, 2, 7, 3, 1]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 3: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 4: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 5: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 6: Crop the matrix to rows 5-6 and columns 1-5 (1-indexed)
matrix = matrix[4:6, :5]

# Step 7: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 8: Horizontally mirror the matrix
matrix = matrix[:, ::-1]

# Step 9: Mirror along the counterdiagonal
matrix = np.rot90(matrix, k=-1).T

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))