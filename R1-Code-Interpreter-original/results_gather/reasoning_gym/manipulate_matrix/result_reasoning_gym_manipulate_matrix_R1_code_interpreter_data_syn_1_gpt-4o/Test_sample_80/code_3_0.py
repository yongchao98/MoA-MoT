import numpy as np

# Initial matrix
matrix = np.array([
    [7, 6, 2, 5, 3, 8],
    [7, 0, 6, 9, 3, 7],
    [7, 4, 5, 3, 3, 9],
    [2, 7, 8, 5, 6, 9],
    [6, 4, 0, 1, 1, 3],
    [4, 6, 1, 9, 0, 1],
    [3, 4, 8, 7, 0, 4],
    [6, 8, 2, 4, 8, 6],
    [2, 8, 9, 1, 3, 9]
])

# Step 1: Identity transformation (no change)
# Step 2: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 3: Set all elements divisible by 3 to zero
matrix[matrix % 3 == 0] = 0

# Step 4: Mirror the matrix along the counterdiagonal
rows, cols = matrix.shape
result = np.zeros((cols, rows), dtype=int)
for i in range(rows):
    for j in range(cols):
        result[j, rows - i - 1] = matrix[i, j]
matrix = result

# Step 5: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 6: Set all elements divisible by 2 to zero
matrix[matrix % 2 == 0] = 0

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))