import numpy as np

# Initial matrix
matrix = np.array([
    [5, 0, 9, 0, 8, 0, 0, 4],
    [5, 0, 1, 7, 7, 1, 5, 3],
    [4, 2, 4, 6, 8, 3, 2, 1],
    [3, 9, 8, 2, 1, 2, 1, 1],
    [8, 8, 4, 9, 1, 0, 7, 1]
])

# Step 2: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 3: Remove every 3rd column (1-indexed)
matrix = np.delete(matrix, np.arange(2, matrix.shape[1], 3), axis=1)

# Step 4: Map each occurrence of 3 to 5
matrix[matrix == 3] = 5

# Step 5: Rotate the matrix 180 degrees
matrix = np.rot90(matrix, 2)

# Print the final matrix
print(matrix)