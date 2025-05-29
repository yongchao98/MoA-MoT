import numpy as np

# Initial matrix
matrix = np.array([
    [3, 2, 7, 7, 3, 2, 2, 7, 1],
    [4, 8, 2, 5, 3, 7, 6, 0, 6]
])

# Step 1: Identity transformation (no change)
# Step 2: Crop to rows 2-2 and columns 9-9 (1-indexed)
matrix = matrix[1:2, 8:9]

# Step 3: Crop to rows 1-1 and columns 1-1 (1-indexed)
matrix = matrix[0:1, 0:1]

# Step 4: Set all elements divisible by 7 to zero
matrix[matrix % 7 == 0] = 0

# Step 5: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 6: Rotate the matrix 360 degrees (no change)
# Step 7: Vertically mirror the matrix (flip upside down)
matrix = np.flipud(matrix)

# Step 8: Map each occurrence of 2 to 3
matrix[matrix == 2] = 3

# Output the final matrix
print(matrix)