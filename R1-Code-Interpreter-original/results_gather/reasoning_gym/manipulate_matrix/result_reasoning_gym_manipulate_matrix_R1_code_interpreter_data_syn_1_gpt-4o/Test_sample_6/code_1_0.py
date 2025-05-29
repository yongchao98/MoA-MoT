import numpy as np

# Original matrix
matrix = np.array([
    [1, 1, 5, 7, 5, 0, 9, 2],
    [9, 5, 5, 4, 5, 9, 4, 7],
    [9, 6, 3, 9, 1, 8, 7, 2],
    [6, 8, 3, 0, 6, 5, 8, 9],
    [8, 0, 2, 7, 1, 4, 1, 2],
    [7, 1, 1, 6, 4, 2, 5, 7],
    [6, 7, 4, 4, 2, 6, 4, 2],
    [5, 1, 3, 8, 9, 3, 4, 0]
])

# Step 2: Remove every 5-th row (1-indexed)
matrix = np.delete(matrix, 4, axis=0)

# Step 3: Remove every 3-th column (1-indexed)
matrix = np.delete(matrix, [2, 5], axis=1)

# Step 4: Crop the matrix to rows 7-7 and columns 1-6 (1-indexed)
final_matrix = matrix[6:7, 0:6]

# Print the final matrix
print(final_matrix)