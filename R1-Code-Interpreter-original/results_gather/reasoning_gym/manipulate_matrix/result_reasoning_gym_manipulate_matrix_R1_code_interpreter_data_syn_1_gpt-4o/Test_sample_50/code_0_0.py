import numpy as np

# Initial matrix
matrix = np.array([
    [8, 4, 2, 2, 4, 2, 0, 6, 7],
    [9, 1, 3, 4, 8, 2, 1, 7, 1],
    [5, 5, 7, 1, 7, 2, 0, 3, 5],
    [4, 0, 8, 5, 8, 7, 6, 6, 4]
])

# Step 2: Remove every 6-th column (1-indexed)
matrix = np.delete(matrix, 5, axis=1)

# Step 3: Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=1)

# Step 4: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 5: Set all elements divisible by 4 to zero
matrix[matrix % 4 == 0] = 0

# Step 6: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 7: Vertically mirror the matrix (cancels out the previous step)
matrix = np.flipud(matrix)

# Step 8: Rotate the matrix 90 degrees (clockwise)
matrix = np.rot90(matrix, k=-1)

# Print the final matrix
print(matrix)