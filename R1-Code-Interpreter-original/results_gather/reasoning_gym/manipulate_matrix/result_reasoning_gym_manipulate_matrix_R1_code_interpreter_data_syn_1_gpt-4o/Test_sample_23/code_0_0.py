import numpy as np

# Initial matrix
matrix = np.array([
    [2, 4, 1, 2, 4],
    [1, 5, 7, 8, 6],
    [6, 8, 6, 5, 6]
])

# Step 2: Crop the matrix to rows 1-3 and columns 3-3 (1-indexed)
matrix = matrix[:, 2:3]

# Step 4: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, 1, axis=0)

# Final output
print(matrix)