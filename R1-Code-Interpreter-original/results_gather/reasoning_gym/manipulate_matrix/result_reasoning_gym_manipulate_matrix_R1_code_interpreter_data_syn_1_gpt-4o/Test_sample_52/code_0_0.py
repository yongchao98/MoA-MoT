import numpy as np

# Initial matrix
matrix = np.array([
    [7, 0],
    [8, 5],
    [3, 7],
    [7, 8],
    [5, 9],
    [2, 8]
])

# Step 1: Identity transformation (no change)
# Step 2: Map each occurrence of 6 to 5 (no 6s in the matrix)

# Step 3: Set all elements divisible by 3 to zero
matrix[matrix % 3 == 0] = 0

# Step 4: Set all elements divisible by 8 to zero
matrix[matrix % 8 == 0] = 0

# Step 5: Rotate the matrix 90 degrees
matrix = np.rot90(matrix, -1)

# Step 6: Horizontally mirror the matrix
matrix = np.fliplr(matrix)

# Print the final matrix
print(matrix)