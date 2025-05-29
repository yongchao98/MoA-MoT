import numpy as np

# Initial matrix
matrix = np.array([
    [8, 8],
    [7, 0],
    [4, 1],
    [3, 7],
    [5, 6],
    [7, 4],
    [3, 7],
    [4, 7],
    [2, 2]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 3: Mirror along the diagonal again (transpose back)
matrix = matrix.T

# Step 4: Vertically mirror the matrix (flip upside down)
matrix = np.flipud(matrix)

# Step 5: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 6: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, np.arange(1, matrix.shape[0], 2), axis=0)

# Prepare the final output in the required format
output = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(f"<<<{output}>>>")