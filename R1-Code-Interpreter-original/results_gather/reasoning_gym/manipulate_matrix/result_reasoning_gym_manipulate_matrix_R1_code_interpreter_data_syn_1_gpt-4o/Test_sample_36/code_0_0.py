import numpy as np

# Initial matrix
matrix = np.array([
    [1, 9, 9, 3, 6, 7, 4, 5, 7, 8],
    [5, 2, 8, 1, 0, 8, 7, 1, 5, 0],
    [1, 9, 2, 3, 3, 0, 8, 6, 7, 8],
    [6, 8, 2, 8, 0, 9, 2, 9, 9, 5],
    [2, 1, 0, 7, 6, 9, 4, 1, 3, 6],
    [6, 6, 6, 1, 4, 6, 2, 0, 0, 6],
    [6, 2, 3, 1, 5, 5, 4, 9, 9, 3],
    [7, 6, 7, 7, 8, 8, 7, 6, 4, 0],
    [0, 8, 1, 7, 3, 9, 4, 9, 3, 5]
])

# Remove every 3rd row
matrix = np.delete(matrix, [2, 5, 8], axis=0)

# Crop to rows 1-5 and columns 5-9
matrix = matrix[:5, 4:9]

# Crop to rows 3-3 and columns 3-3
matrix = matrix[2:3, 2:3]

# Map each occurrence of 6 to 8
matrix[matrix == 6] = 8

# Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Mirror the matrix along the diagonal (transpose again)
matrix = matrix.T

# Crop to rows 1-1 and columns 1-1
matrix = matrix[0:1, 0:1]

# Vertically mirror the matrix
matrix = np.flipud(matrix)

# Output the final matrix
print(matrix)