import numpy as np

# Initial matrix
matrix = np.array([
    [9, 4, 1, 5, 3],
    [7, 9, 2, 1, 4],
    [7, 3, 8, 6, 6],
    [2, 5, 4, 7, 1],
    [6, 1, 8, 7, 9],
    [1, 8, 9, 7, 4],
    [5, 9, 0, 8, 5],
    [6, 6, 2, 5, 6],
    [1, 7, 1, 9, 0],
    [0, 3, 7, 6, 4]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 3: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 4: Set all elements divisible by 6 to zero
matrix[matrix % 6 == 0] = 0

# Step 5: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 6: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 7: Mirror along the diagonal
matrix = matrix.T

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))