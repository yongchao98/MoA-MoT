import numpy as np

# Initial matrix
matrix = np.array([
    [1, 5, 4, 0, 4, 5],
    [0, 9, 7, 8, 9, 5],
    [9, 8, 2, 9, 8, 0],
    [3, 9, 2, 6, 4, 1]
])

# Step 2: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, [1, 3], axis=0)

# Step 3: Map each occurrence of 5 to 9
matrix[matrix == 5] = 9

# Step 4: Crop the matrix to rows 2-2 and columns 1-4 (1-indexed)
matrix = matrix[1:2, 0:4]

# Step 5: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 6: Set all elements divisible by 1 to zero
matrix = np.zeros_like(matrix)

# Step 7: Remove every 3-th column (1-indexed)
matrix = np.delete(matrix, [2], axis=1)

# Step 8: Rotate the matrix 90 degrees
matrix = np.rot90(matrix, k=1)

# Step 9: Rotate the matrix 180 degrees
matrix = np.rot90(matrix, k=2)

# Step 10: Mirror the matrix along the diagonal
matrix = matrix.T

# Step 11: Crop the matrix to rows 1-1 and columns 1-2 (1-indexed)
matrix = matrix[0:1, 0:2]

# Print the final matrix
print(matrix)