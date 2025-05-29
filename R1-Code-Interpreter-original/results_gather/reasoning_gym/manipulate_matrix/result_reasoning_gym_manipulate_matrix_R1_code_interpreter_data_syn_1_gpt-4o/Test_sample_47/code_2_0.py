import numpy as np

# Initial matrix
matrix = np.array([
    [5, 9, 2, 3, 3, 7, 2, 1, 5],
    [0, 3, 1, 7, 4, 2, 6, 7, 5],
    [9, 1, 6, 1, 4, 7, 0, 5, 7],
    [1, 0, 7, 1, 1, 3, 2, 6, 4]
])

# Step 1: Identity transformation (no change)
# Step 2: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 3: Set all elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Step 4: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 5: Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=3)

# Step 6: Remove every 2-nd column (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[1], 2), axis=1)

# Step 7: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Print the final matrix
print(matrix)