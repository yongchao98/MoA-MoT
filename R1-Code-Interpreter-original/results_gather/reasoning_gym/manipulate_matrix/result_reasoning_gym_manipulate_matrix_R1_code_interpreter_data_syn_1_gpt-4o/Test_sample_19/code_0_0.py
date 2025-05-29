import numpy as np

# Initial matrix
matrix = np.array([
    [8, 2, 3, 3, 8, 4, 6, 8, 7, 8],
    [8, 5, 9, 5, 6, 4, 6, 0, 1, 0],
    [8, 4, 0, 7, 2, 5, 9, 4, 5, 7],
    [3, 0, 7, 2, 7, 0, 0, 1, 9, 1],
    [2, 8, 7, 0, 6, 6, 3, 5, 1, 7]
])

# Step 2: Crop the matrix to rows 5-5 and columns 8-9 (1-indexed)
matrix = matrix[4:5, 7:9]

# Step 3: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 4: Crop the matrix to rows 1-1 and columns 1-2 (1-indexed)
matrix = matrix[0:1, 0:2]

# Step 5: Remove every 2-nd column (1-indexed)
matrix = matrix[:, ::2]

# Step 6: Map each occurrence of 4 to 8
matrix[matrix == 4] = 8

# Step 7: Rotate the matrix 360 degrees (no change)
# No operation needed

# Step 8: Map each occurrence of 9 to 6
matrix[matrix == 9] = 6

# Step 9: Map each occurrence of 4 to 3
matrix[matrix == 4] = 3

# Step 10: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Print the final matrix
print(matrix)