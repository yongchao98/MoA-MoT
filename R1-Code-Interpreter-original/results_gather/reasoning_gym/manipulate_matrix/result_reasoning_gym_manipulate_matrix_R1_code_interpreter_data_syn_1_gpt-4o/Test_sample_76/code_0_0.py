import numpy as np

# Initial matrix
matrix = np.array([
    [6, 4, 6],
    [3, 1, 4],
    [3, 2, 3],
    [5, 3, 1]
])

# Step 1: Identity transformation (no change)
# Step 2: Rotate 360 degrees (no change)

# Step 3: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 4: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 5: Remove every 3rd column (1-indexed)
matrix = np.delete(matrix, 2, axis=1)

# Print the final matrix
print(matrix)