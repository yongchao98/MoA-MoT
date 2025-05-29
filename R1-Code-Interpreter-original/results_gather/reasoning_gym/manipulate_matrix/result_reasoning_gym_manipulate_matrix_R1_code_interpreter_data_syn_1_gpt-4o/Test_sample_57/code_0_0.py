import numpy as np

# Initial matrix
matrix = np.array([
    [2, 5, 1, 2, 6, 5],
    [3, 4, 5, 1, 9, 2],
    [7, 9, 0, 3, 2, 9]
])

# Step 2: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, 1, axis=0)

# Step 3: Map each occurrence of 1 to 9
matrix[matrix == 1] = 9

# Step 4: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(matrix).T

# Step 5: Rotate the matrix 270 degrees
matrix = np.rot90(matrix, k=3)

# Step 6: Crop the matrix to rows 1-1 and columns 1-2 (1-indexed)
matrix = matrix[0:1, 0:2]

# Step 7: Map each occurrence of 1 to 2 (no 1s present)

# Step 8: Mirror the matrix along the diagonal
matrix = matrix.T

# Step 9: Crop the matrix to rows 2-2 and columns 1-1 (1-indexed)
matrix = matrix[1:2, 0:1]

# Output the final matrix
print(matrix)