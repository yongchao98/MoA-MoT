import numpy as np

# Initial matrix
matrix = np.array([
    [8, 1, 2, 6, 3, 4, 0, 3, 1, 9],
    [0, 1, 2, 8, 4, 6, 9, 6, 5, 5],
    [1, 5, 4, 9, 2, 1, 8, 1, 9, 1],
    [4, 5, 1, 4, 0, 5, 6, 1, 7, 7],
    [3, 3, 2, 4, 3, 0, 0, 6, 0, 5],
    [5, 7, 7, 9, 8, 2, 3, 7, 7, 5],
    [9, 0, 4, 2, 0, 3, 0, 9, 9, 8],
    [8, 4, 5, 9, 3, 6, 1, 5, 5, 1]
])

# Step 2: Set all elements divisible by 1 to zero
matrix[:] = 0

# Step 3: Remove every 2-nd column (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[1], 2), axis=1)

# Step 4: Set all elements divisible by 8 to zero (no effect since all are zero)

# Step 5: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Print the final matrix
print(matrix)