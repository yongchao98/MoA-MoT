import numpy as np

# Initial matrix
matrix = np.array([
    [9, 4, 9, 2, 5],
    [2, 9, 4, 9, 0],
    [1, 2, 7, 8, 6],
    [1, 4, 4, 2, 7],
    [2, 7, 2, 0, 4]
])

# Step 2: Set all elements divisible by 7 to zero
matrix[matrix % 7 == 0] = 0

# Step 3: Rotate the matrix 180 degrees
matrix = np.rot90(matrix, 2)

# Step 4: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[0], 2), axis=0)

# Step 5: Rotate the matrix 90 degrees
matrix = np.rot90(matrix, -1)

# Step 6: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 8: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix))

# Step 9: Set all elements divisible by 6 to zero
matrix[matrix % 6 == 0] = 0

# Step 10: Map each occurrence of 6 to 3
matrix[matrix == 6] = 3

# Step 11: Remove every 5-th row (1-indexed)
matrix = np.delete(matrix, np.arange(4, matrix.shape[0], 5), axis=0)

# Print the final matrix
print(matrix)