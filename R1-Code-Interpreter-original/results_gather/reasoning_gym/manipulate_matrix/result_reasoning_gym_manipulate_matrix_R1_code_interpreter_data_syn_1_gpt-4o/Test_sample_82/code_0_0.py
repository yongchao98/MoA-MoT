import numpy as np

# Initial matrix
matrix = np.array([
    [5, 3, 9],
    [0, 2, 8],
    [8, 8, 4],
    [9, 1, 8],
    [5, 5, 1]
])

# Step 1: Identity transformation (no change)
# Step 2: Remove every 5-th row (1-indexed)
matrix = np.delete(matrix, 4, axis=0)

# Step 3: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Print the final matrix
print(matrix)