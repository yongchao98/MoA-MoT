import numpy as np

# Original matrix
matrix = np.array([
    [3, 0, 3, 7, 6],
    [5, 4, 4, 7, 9],
    [6, 4, 9, 1, 0],
    [7, 3, 7, 1, 7],
    [8, 6, 5, 0, 1]
])

# Identity transformation (no change)
identity_matrix = matrix.copy()

# Mirror the matrix along the diagonal (transpose)
transposed_matrix = identity_matrix.T

# Horizontally mirror the matrix (reverse the order of rows)
horizontally_mirrored_matrix = transposed_matrix[::-1]

# Vertically mirror the matrix (reverse the order of elements in each row)
final_matrix = np.array([row[::-1] for row in horizontally_mirrored_matrix])

# Print the final matrix
for row in final_matrix:
    print(' '.join(map(str, row)))