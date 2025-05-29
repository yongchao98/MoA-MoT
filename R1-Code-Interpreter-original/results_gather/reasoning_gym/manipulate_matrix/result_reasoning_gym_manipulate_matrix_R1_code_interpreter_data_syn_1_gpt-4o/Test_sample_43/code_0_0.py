import numpy as np

# Initial matrix
matrix = np.array([
    [7, 0, 8, 7, 5, 1, 5, 1],
    [4, 4, 7, 5, 9, 7, 4, 6],
    [0, 6, 8, 9, 7, 3, 2, 3],
    [2, 4, 7, 9, 2, 7, 0, 4],
    [9, 0, 1, 0, 4, 5, 3, 0],
    [1, 6, 1, 4, 5, 8, 1, 7],
    [4, 4, 9, 3, 3, 6, 0, 0],
    [4, 5, 7, 4, 0, 8, 2, 8],
    [8, 7, 0, 4, 8, 5, 2, 5]
])

# Step 2: Remove every 9-th row (1-indexed)
matrix = np.delete(matrix, 8, axis=0)

# Step 3: Crop the matrix to rows 2-5 and columns 5-7 (1-indexed)
matrix = matrix[1:5, 4:7]

# Step 4: Remove every 3-th column (1-indexed)
matrix = np.delete(matrix, 2, axis=1)

# Step 5: Remove every 3-th row (1-indexed)
matrix = np.delete(matrix, 2, axis=0)

# Step 6: Rotate the matrix 360 degrees (no change)
# Final matrix
print(matrix)