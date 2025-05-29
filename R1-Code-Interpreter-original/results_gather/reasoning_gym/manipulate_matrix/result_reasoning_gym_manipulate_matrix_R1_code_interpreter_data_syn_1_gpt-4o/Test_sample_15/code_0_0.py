import numpy as np

# Initial matrix
matrix = np.array([
    [6, 8, 7],
    [9, 5, 4],
    [8, 6, 6],
    [9, 0, 6],
    [0, 1, 0],
    [4, 7, 8],
    [0, 5, 9]
])

# Step 2: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 3: Map each occurrence of 5 to 7
matrix[matrix == 5] = 7

# Step 4: Remove every 2-nd column (1-indexed)
matrix = np.delete(matrix, 1, axis=1)

# Step 5: Map each occurrence of 5 to 9 (no effect as there are no 5s)
# Step 6: Map each occurrence of 3 to 9 (no effect as there are no 3s)

# Step 7: Crop the matrix to rows 6-6 and columns 1-2 (1-indexed)
matrix = matrix[5:6, 0:2]

# Step 8: Rotate the matrix 360 degrees (no change)
# Final matrix
print(matrix)