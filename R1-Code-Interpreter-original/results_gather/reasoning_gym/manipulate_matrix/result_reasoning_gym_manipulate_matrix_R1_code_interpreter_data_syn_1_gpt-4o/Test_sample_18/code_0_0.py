import numpy as np

# Initial matrix
matrix = np.array([
    [5, 3, 7, 5, 1, 8],
    [5, 3, 8, 1, 2, 5],
    [1, 5, 5, 3, 3, 1],
    [5, 9, 9, 6, 1, 1],
    [9, 5, 6, 0, 3, 5],
    [4, 6, 1, 3, 7, 4],
    [9, 7, 2, 4, 9, 3],
    [8, 9, 1, 6, 9, 1],
    [8, 7, 6, 8, 8, 2]
])

# Step 2: Map each occurrence of 8 to 4
matrix[matrix == 8] = 4

# Step 3: Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=1)

# Step 4: Crop the matrix to rows 2-2 and columns 6-6 (1-indexed)
matrix = matrix[1:2, 5:6]

# Step 5: Set all elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Step 6: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 7: Crop the matrix to rows 1-1 and columns 1-1 (1-indexed)
matrix = matrix[0:1, 0:1]

# Step 8: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 9: Set all elements divisible by 1 to zero
matrix[matrix % 1 == 0] = 0

# Output the final matrix
print(matrix)