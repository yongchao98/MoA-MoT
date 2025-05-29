import numpy as np

# Initial matrix
matrix = np.array([
    [4, 7, 8, 8, 2],
    [9, 8, 5, 8, 7],
    [7, 4, 1, 4, 6],
    [4, 8, 7, 0, 9]
])

# Step 2: Remove every 3rd row
matrix = np.delete(matrix, 2, axis=0)

# Step 3: Remove every 3rd column
matrix = np.delete(matrix, 2, axis=1)

# Step 4: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Step 5: Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=1)

# Step 6: Remove every 2nd column
matrix = np.delete(matrix, 1, axis=1)

# Step 7: Crop the matrix to rows 1-3 and columns 1-2
matrix = matrix[:3, :2]

# Step 8: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Print the final matrix
print(matrix)