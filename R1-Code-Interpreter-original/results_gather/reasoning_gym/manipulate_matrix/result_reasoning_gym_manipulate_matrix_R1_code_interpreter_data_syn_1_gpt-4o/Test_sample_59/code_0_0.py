import numpy as np

# Initial matrix
matrix = np.array([
    [2, 6, 7, 7],
    [7, 7, 4, 3],
    [3, 9, 5, 7]
])

# Step 1: Identity transformation (no change)
# matrix remains the same

# Step 2: Remove every 2-nd column (1-indexed)
matrix = np.delete(matrix, [1, 3], axis=1)

# Step 3: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 4: Rotate the matrix 90 degrees clockwise
matrix = np.rot90(matrix, -1)

# Step 5: Remove every 3-th row (1-indexed)
matrix = np.delete(matrix, 2, axis=0)

# Print the final matrix
print(matrix)