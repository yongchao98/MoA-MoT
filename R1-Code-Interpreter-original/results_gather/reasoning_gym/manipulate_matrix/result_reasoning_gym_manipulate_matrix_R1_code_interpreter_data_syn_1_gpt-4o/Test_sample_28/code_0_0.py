import numpy as np

# Initial matrix
matrix = np.array([
    [9, 2, 4, 3, 7, 7],
    [5, 0, 1, 2, 2, 6],
    [0, 3, 0, 6, 2, 2],
    [4, 9, 2, 8, 9, 8],
    [9, 0, 9, 9, 1, 4],
    [2, 0, 9, 3, 4, 3]
])

# Step 2: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 3: Remove every 3rd column (1-indexed)
matrix = np.delete(matrix, [2, 5], axis=1)

# Step 4: Crop the matrix to rows 3-4 and columns 2-4 (1-indexed)
matrix = matrix[2:4, 1:4]

# Step 5: Set all elements divisible by 2 to zero
matrix[matrix % 2 == 0] = 0

# Step 6: Map each occurrence of 7 to 5
matrix[matrix == 7] = 5

# Step 7: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 8: Remove every 2nd row (1-indexed)
matrix = np.delete(matrix, [1], axis=0)

# Step 9: Rotate the matrix 90 degrees clockwise
matrix = np.rot90(matrix, -1)

# Step 10: Rotate the matrix 360 degrees (no change)
# No operation needed

# Print the final matrix
print(matrix)