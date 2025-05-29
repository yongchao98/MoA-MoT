import numpy as np

# Initial matrix
matrix = np.array([
    [9, 8, 8, 9, 4, 1, 9, 8, 7, 2],
    [1, 9, 5, 0, 0, 0, 2, 9, 2, 9],
    [9, 2, 6, 7, 2, 1, 2, 7, 1, 6],
    [5, 6, 0, 4, 2, 2, 3, 9, 2, 3],
    [9, 2, 7, 0, 6, 4, 4, 4, 1, 9],
    [8, 0, 5, 6, 2, 9, 9, 4, 5, 1],
    [1, 7, 7, 2, 4, 5, 2, 5, 3, 8],
    [8, 7, 3, 2, 5, 0, 0, 1, 2, 8]
])

# Step 2: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 3: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix))

# Step 4: Rotate 90 degrees (transpose and reverse each row)
matrix = np.rot90(matrix, k=-1)

# Step 5: Remove every 10-th row (no effect here)

# Step 6: Mirror along the diagonal (transpose)
matrix = matrix.T

# Step 7: Map each occurrence of 1 to 0
matrix[matrix == 1] = 0

# Step 8: Set all elements divisible by 2 to zero
matrix[matrix % 2 == 0] = 0

# Step 9: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix))

# Step 10: Crop to rows 4-9 and columns 4-6 (1-indexed)
matrix = matrix[3:9, 3:6]

# Step 11: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix))

# Output the final matrix
output = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(f"<<<{output}>>>")