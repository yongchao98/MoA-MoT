import numpy as np

# Initial matrix
matrix = np.array([
    [0, 7, 8, 5, 3, 6],
    [0, 1, 4, 2, 9, 4]
])

# Step 2: Remove every 4-th column (1-indexed)
matrix = np.delete(matrix, 3, axis=1)

# Step 3: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 4: Mirror the matrix along the diagonal
matrix = matrix.T

# Step 5: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 6: Remove every 2-nd column (1-indexed)
matrix = np.delete(matrix, 1, axis=1)

# Step 7: Map each occurrence of 2 to 3
matrix[matrix == 2] = 3

# Step 8: Set all elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Step 9: Set all elements divisible by 3 to zero
matrix[matrix % 3 == 0] = 0

# Print the final matrix
print(matrix)