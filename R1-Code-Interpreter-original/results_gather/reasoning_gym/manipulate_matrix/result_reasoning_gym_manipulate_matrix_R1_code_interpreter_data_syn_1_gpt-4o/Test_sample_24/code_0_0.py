import numpy as np

# Original matrix
matrix = np.array([
    [7, 8],
    [4, 8],
    [6, 3],
    [1, 5],
    [9, 1],
    [5, 4],
    [9, 3],
    [0, 9],
    [3, 0],
    [9, 8]
])

# Step 1: Identity transformation (no change)
identity_matrix = matrix

# Step 2: Mirror the matrix along the diagonal (transpose)
transposed_matrix = identity_matrix.T

# Step 3: Rotate the matrix 180 degrees
rotated_matrix = np.rot90(transposed_matrix, 2)

# Step 4: Mirror the matrix along the diagonal again (transpose)
final_matrix = rotated_matrix.T

# Print the final matrix
for row in final_matrix:
    print(' '.join(map(str, row)))