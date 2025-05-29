import numpy as np

# Initial matrix
matrix = np.array([
    [4, 4, 7, 1, 6, 1, 1, 3, 4],
    [1, 4, 9, 0, 3, 2, 5, 3, 4],
    [9, 0, 0, 9, 2, 3, 9, 8, 6],
    [6, 5, 7, 3, 1, 0, 0, 5, 9],
    [5, 1, 0, 9, 2, 9, 0, 6, 7],
    [0, 5, 0, 5, 0, 5, 6, 1, 4],
    [8, 8, 1, 3, 8, 3, 9, 6, 7],
    [2, 0, 5, 3, 6, 5, 8, 2, 4],
    [7, 1, 6, 2, 5, 2, 2, 4, 0]
])

# Step 2: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 3: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 4: Set all elements divisible by 7 to zero
matrix[matrix % 7 == 0] = 0

# Step 5: Remove every 6-th row (1-indexed)
matrix = np.delete(matrix, np.arange(5, matrix.shape[0], 6), axis=0)

# Step 6: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 7: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[0], 2), axis=0)

# Step 8: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 9: Crop the matrix to rows 1-4 and columns 8-8 (1-indexed)
final_matrix = matrix[:4, 7:8]

# Output the final matrix
print(final_matrix)