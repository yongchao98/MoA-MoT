import numpy as np

# Initial matrix
matrix = np.array([
    [0, 4],
    [1, 0],
    [8, 2],
    [3, 4],
    [6, 5],
    [8, 4],
    [8, 1]
])

# Step 1: Identity transformation (no change)
# Step 2: Mirror along the counterdiagonal (swap rows)
matrix = matrix[::-1]

# Step 3: Mirror along the counterdiagonal (swap rows again)
matrix = matrix[::-1]

# Step 4: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 5: Crop to rows 2-2 and columns 4-5 (1-indexed)
# Adjusted for transposed matrix, this should be columns 1-1 and rows 4-5 (0-indexed)
matrix = matrix[3:5, 1:2]

# Step 6: Horizontally mirror the matrix (no effect on empty matrix)
# Step 7: Mirror along the counterdiagonal (no effect on empty matrix)
# Step 8: Crop to rows 2-2 and columns 1-1 (1-indexed) (no effect on empty matrix)
# Step 9: Crop to rows 1-1 and columns 1-1 (1-indexed) (no effect on empty matrix)
# Step 10: Horizontally mirror the matrix (no effect on empty matrix)
# Step 11: Mirror along the counterdiagonal (no effect on empty matrix)

# Output the final matrix
print(matrix)