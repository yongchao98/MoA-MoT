import numpy as np

# Initial matrix
matrix = np.array([
    [1, 4, 4, 6],
    [0, 7, 5, 4],
    [6, 9, 4, 8],
    [0, 6, 0, 6],
    [5, 7, 5, 3],
    [0, 6, 7, 2],
    [4, 2, 1, 4],
    [7, 4, 1, 5],
    [6, 5, 0, 8],
    [0, 4, 3, 0]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 3: Vertically mirror the matrix (reverse rows)
matrix = matrix[::-1]

# Step 4: Rotate the matrix 180 degrees (reverse both rows and columns)
matrix = matrix[::-1, ::-1]

# Step 5: Set all elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Step 6: Map each occurrence of 3 to 7
matrix[matrix == 3] = 7

# Step 7: Rotate the matrix 180 degrees (reverse both rows and columns)
matrix = matrix[::-1, ::-1]

# Step 8: Map each occurrence of 8 to 1
matrix[matrix == 8] = 1

# Step 9: Set all elements divisible by 5 to zero
matrix[matrix % 5 == 0] = 0

# Print the final matrix
print(matrix)