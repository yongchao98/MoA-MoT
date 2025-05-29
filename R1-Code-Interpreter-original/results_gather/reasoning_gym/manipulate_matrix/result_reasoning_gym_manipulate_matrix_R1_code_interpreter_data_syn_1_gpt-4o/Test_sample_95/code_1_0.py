import numpy as np

# Initial matrix
matrix = np.array([[0, 0, 4],
                   [3, 6, 3]])

# Step 1: Identity transformation (no change)
# Step 2: Set all elements divisible by 5 to zero
matrix[matrix % 5 == 0] = 0

# Step 3: Mirror the matrix along the diagonal
matrix = matrix.T

# Step 4: Rotate the matrix 270 degrees (90 degrees counterclockwise)
matrix = np.rot90(matrix, k=1)

# Step 5: Remove every 2-nd column (1-indexed)
matrix = np.delete(matrix, 1, axis=1)

# Step 6: Remove every 2-nd column (1-indexed)
# Since there's only one column left, this step doesn't change the matrix

# Print the final matrix
print(matrix)