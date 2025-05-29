import numpy as np

# Initial matrix
matrix = np.array([
    [0, 6, 3, 5, 3, 1],
    [8, 2, 2, 4, 9, 3],
    [2, 7, 9, 7, 2, 1],
    [4, 1, 4, 1, 6, 0],
    [3, 9, 2, 1, 3, 7],
    [0, 4, 9, 4, 5, 2],
    [5, 7, 8, 9, 1, 8],
    [2, 3, 2, 8, 7, 2],
    [9, 5, 4, 9, 5, 9],
    [7, 5, 5, 2, 8, 2]
])

# Step 2: Remove every 5-th row
matrix = np.delete(matrix, [4, 9], axis=0)

# Step 3: Rotate the matrix 180 degrees
matrix = np.rot90(matrix, 2)

# Step 4: Vertically mirror the matrix
matrix = np.fliplr(matrix)

# Step 5: Remove every 3-th column
matrix = np.delete(matrix, [2, 5], axis=1)

# Step 6: Map each occurrence of 1 to 6
matrix[matrix == 1] = 6

# Step 7: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 8: Horizontally mirror the matrix again (no change)
matrix = np.fliplr(matrix)

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))