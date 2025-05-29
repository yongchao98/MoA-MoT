import numpy as np

# Original matrix
matrix = np.array([
    [6, 1, 5, 8, 8, 0, 0],
    [4, 2, 7, 0, 8, 8, 3],
    [9, 4, 9, 5, 6, 0, 3]
])

# Identity transformation (no change)
identity_matrix = matrix.copy()

# Rotate 360 degrees (no change)
rotated_matrix = identity_matrix.copy()

# Mirror along the diagonal (transpose)
mirrored_once = rotated_matrix.T

# Mirror along the diagonal again (transpose again)
final_matrix = mirrored_once.T

# Print the final matrix
print(final_matrix)