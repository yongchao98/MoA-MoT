import numpy as np

# Initial matrix
matrix = np.array([
    [1, 6, 4, 1, 0, 6],
    [8, 8, 5, 5, 0, 2]
])

# Identity transformation (no change)
identity_matrix = matrix.copy()

# Mirror along the diagonal (transpose)
transposed_matrix = identity_matrix.T

# Mirror along the diagonal again (transpose back)
original_matrix = transposed_matrix.T

# Vertically mirror the matrix (flip upside down)
final_matrix = np.flipud(original_matrix)

# Print the final matrix
for row in final_matrix:
    print(' '.join(map(str, row)))