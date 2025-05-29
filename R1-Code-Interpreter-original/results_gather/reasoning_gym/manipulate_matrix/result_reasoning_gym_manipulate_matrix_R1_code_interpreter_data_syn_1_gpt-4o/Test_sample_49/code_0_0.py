import numpy as np

# Initial matrix
matrix = np.array([
    [5, 8, 8, 6, 7, 0, 2, 2, 7],
    [9, 0, 5, 4, 5, 3, 4, 2, 9],
    [1, 3, 8, 4, 9, 8, 2, 5, 2],
    [6, 3, 7, 3, 7, 0, 8, 6, 5],
    [4, 5, 4, 9, 5, 1, 2, 0, 6],
    [4, 4, 2, 3, 7, 1, 3, 3, 1],
    [2, 9, 4, 6, 0, 0, 6, 8, 2],
    [2, 1, 6, 6, 3, 0, 9, 5, 0],
    [3, 5, 7, 6, 4, 1, 6, 0, 6],
    [3, 1, 0, 8, 7, 3, 5, 4, 9]
])

# Step 2: Crop to rows 8-10 and columns 7-9
matrix = matrix[7:10, 6:9]

# Step 3: Remove every 2-nd row
matrix = matrix[::2]

# Step 4: Vertically mirror the matrix
matrix = matrix[::-1]

# Step 5: Remove every 3-th column
matrix = np.delete(matrix, np.s_[2::3], axis=1)

# Step 6: Remove every 2-nd column
matrix = np.delete(matrix, np.s_[1::2], axis=1)

# Step 7: Set all elements divisible by 9 to zero
matrix[matrix % 9 == 0] = 0

# Step 8: Crop to rows 1-1 and columns 1-1
matrix = matrix[0:1, 0:1]

# Step 9: Crop to rows 1-1 and columns 1-1 (redundant)
# Already a single element

# Step 10: Map each occurrence of 2 to 7
matrix[matrix == 2] = 7

# Step 11: Horizontally mirror the matrix (no effect on single element)

# Output the final matrix
print(matrix)