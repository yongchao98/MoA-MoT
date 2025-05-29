import numpy as np

# Initial matrix
matrix = np.array([
    [8, 8, 6, 1],
    [1, 6, 3, 2],
    [4, 4, 3, 7],
    [9, 0, 6, 4],
    [6, 0, 5, 3],
    [2, 7, 0, 8]
])

# Step 2: Crop the matrix to rows 4-4 and columns 1-2 (1-indexed)
matrix = matrix[3:4, 0:2]

# Step 3: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 4: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(matrix).T

# Step 5: Mirror the matrix along the counterdiagonal (revert to original)
matrix = np.fliplr(matrix).T

# Step 6: Remove every 2-nd column (1-indexed)
matrix = matrix[:, ::2]

# Step 7: Crop the matrix to rows 1-1 and columns 1-1 (1-indexed)
matrix = matrix[0:1, 0:1]

# Step 8: Set all elements divisible by 6 to zero
matrix[matrix % 6 == 0] = 0

# Step 9: Crop the matrix to rows 1-1 and columns 1-1 (1-indexed) (redundant)

# Step 10: Map each occurrence of 3 to 9
matrix[matrix == 3] = 9

# Output the final matrix
print(matrix)