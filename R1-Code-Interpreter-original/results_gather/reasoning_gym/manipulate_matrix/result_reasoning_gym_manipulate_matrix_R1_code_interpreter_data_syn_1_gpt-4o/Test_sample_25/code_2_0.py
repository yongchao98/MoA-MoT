import numpy as np

# Initial matrix
matrix = np.array([
    [9, 9, 0, 3, 8, 5, 7],
    [4, 6, 6, 0, 2, 8, 9]
])

# Step 1: Identity transformation (no change)
# Step 2: Remove every 2-nd column (1-indexed)
matrix = np.delete(matrix, [1, 3, 5], axis=1)

# Step 3: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 4: Map each occurrence of 9 to 2
matrix[matrix == 9] = 2

# Step 5: Map each occurrence of 2 to 6
matrix[matrix == 2] = 6

# Format the output
output = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(f"<<<{output}>>>")