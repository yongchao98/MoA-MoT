import numpy as np

# Initial matrix
matrix = np.array([
    [8, 1, 6, 4, 0, 6, 1],
    [5, 9, 9, 1, 2, 3, 8],
    [7, 7, 3, 3, 8, 0, 2],
    [7, 9, 3, 6, 1, 0, 7],
    [4, 6, 1, 2, 1, 2, 1],
    [9, 5, 1, 8, 9, 6, 3],
    [2, 9, 7, 5, 0, 2, 7]
])

# Step 3: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 4: Crop the matrix to rows 1-4 and columns 5-7 (1-indexed)
cropped_matrix = matrix[0:4, 4:7]

# Step 5: Horizontally mirror the cropped matrix
cropped_matrix = np.fliplr(cropped_matrix)

# Step 6: Mirror the matrix along the counterdiagonal
m, n = cropped_matrix.shape
counterdiagonal_mirrored = np.zeros_like(cropped_matrix)
for i in range(m):
    for j in range(n):
        counterdiagonal_mirrored[j, i] = cropped_matrix[i, j]

# Print the final matrix
for row in counterdiagonal_mirrored:
    print(' '.join(map(str, row)))