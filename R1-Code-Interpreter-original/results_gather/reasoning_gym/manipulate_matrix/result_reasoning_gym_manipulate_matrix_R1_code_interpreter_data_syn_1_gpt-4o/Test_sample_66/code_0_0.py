import numpy as np

# Initial matrix
matrix = np.array([
    [5, 6],
    [0, 7],
    [4, 8],
    [1, 7],
    [3, 6],
    [0, 5],
    [0, 0]
])

# Step 1: Identity transformation (no change)
# Step 2: Crop the matrix to rows 3-3 and columns 2-2 (1-indexed)
matrix = matrix[2:3, 1:2]

# Step 3: Map each occurrence of 6 to 2
matrix[matrix == 6] = 2

# Step 4: Map each occurrence of 3 to 4
matrix[matrix == 3] = 4

# Step 5: Set all elements divisible by 3 to zero
matrix[matrix % 3 == 0] = 0

# Step 6: Vertically mirror the matrix
matrix = np.flipud(matrix)

# Step 7: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 8: Crop the matrix to rows 1-1 and columns 1-1 (1-indexed)
matrix = matrix[0:1, 0:1]

# Output the final matrix
print(matrix)