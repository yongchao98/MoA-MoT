import numpy as np

# Initial matrix
matrix = np.array([
    [3, 2, 6, 9, 4, 8, 4],
    [6, 9, 0, 3, 0, 7, 2],
    [2, 5, 2, 4, 9, 9, 0],
    [1, 9, 0, 6, 4, 3, 1],
    [5, 7, 2, 8, 0, 1, 1],
    [2, 2, 0, 4, 5, 3, 9],
    [0, 9, 9, 2, 6, 0, 4],
    [6, 2, 1, 1, 6, 0, 7]
])

# Step 2: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 3: Rotate the matrix 90 degrees clockwise
matrix = np.rot90(matrix, -1)

# Step 4: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 5: Crop the matrix to rows 7-8 and columns 5-5 (1-indexed)
# Adjusting for 0-indexing, this means rows 6-7 and column 4
matrix = matrix[6:8, 4:5]

# Step 6: Remove every 2-nd row (1-indexed)
matrix = matrix[::2]

# Step 7: Mirror the matrix along the diagonal
matrix = matrix.T

# Step 8: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Print the final matrix
print(matrix)