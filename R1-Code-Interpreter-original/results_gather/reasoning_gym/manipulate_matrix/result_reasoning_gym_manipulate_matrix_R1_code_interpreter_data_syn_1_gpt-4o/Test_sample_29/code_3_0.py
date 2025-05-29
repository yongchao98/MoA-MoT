import numpy as np

# Initial matrix
matrix = np.array([
    [8, 3, 5, 4],
    [5, 0, 7, 5],
    [5, 1, 6, 7],
    [5, 2, 5, 0],
    [8, 1, 2, 5],
    [0, 6, 2, 2],
    [6, 4, 4, 5],
    [7, 5, 9, 5],
    [2, 9, 4, 6]
])

# Step 1: Identity transformation (no change)
# Step 2: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 3: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 4: Map each occurrence of 7 to 0
matrix[matrix == 7] = 0

# Step 5: Map each occurrence of 1 to 4
matrix[matrix == 1] = 4

# Step 6: Remove every 3rd column (1-indexed)
matrix = np.delete(matrix, np.s_[2::3], axis=1)

# Step 7: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 8: Mirror the matrix along the counterdiagonal
matrix = np.rot90(matrix, 2).T

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))