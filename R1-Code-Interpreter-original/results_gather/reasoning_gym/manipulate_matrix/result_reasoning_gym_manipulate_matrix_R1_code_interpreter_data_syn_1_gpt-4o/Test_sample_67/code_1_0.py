import numpy as np

# Initial matrix
matrix = np.array([
    [3, 9, 2, 3, 3, 2, 1, 4],
    [6, 8, 8, 1, 9, 7, 3, 6],
    [9, 9, 3, 5, 6, 3, 2, 7],
    [2, 9, 5, 1, 8, 1, 7, 5],
    [6, 0, 4, 8, 0, 6, 3, 5],
    [8, 1, 7, 8, 3, 9, 7, 0],
    [5, 2, 3, 2, 7, 9, 2, 8],
    [2, 0, 1, 0, 2, 7, 6, 1]
])

# Step 2: Set all elements divisible by 3 to zero
matrix[matrix % 3 == 0] = 0

# Step 3: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 4: Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=1)

# Step 5: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 6: Mirror the matrix along the counterdiagonal
matrix = np.flipud(np.fliplr(matrix))

# Step 7: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[0], 2), axis=0)

# Step 8: Set all elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Prepare the final output
output = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(f"<<<{output}>>>")