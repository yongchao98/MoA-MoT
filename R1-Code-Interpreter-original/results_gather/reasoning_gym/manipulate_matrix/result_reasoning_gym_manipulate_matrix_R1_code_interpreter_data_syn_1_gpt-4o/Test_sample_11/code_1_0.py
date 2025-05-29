import numpy as np

# Initial matrix
matrix = np.array([
    [2, 7, 5, 1, 7, 9],
    [7, 9, 0, 8, 6, 9]
])

# Step 2: Set all elements divisible by 9 to zero
matrix[matrix % 9 == 0] = 0

# Step 3: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, 1, axis=0)

# Step 4: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 5: Rotate the matrix 90 degrees clockwise
matrix = np.rot90(matrix, -1)

# Step 6: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 7: Rotate the matrix 360 degrees (no change)
# No operation needed as 360-degree rotation results in the same matrix

# Print the final matrix
print(matrix)