import numpy as np

# Initial matrix
matrix = np.array([[2, 1, 9],
                   [5, 5, 8]])

# Step 1: Identity transformation (no change)
# Step 2, 3, 4: Mirror along the counterdiagonal three times
# Mirroring along the counterdiagonal is equivalent to reversing the order of rows and columns
matrix = np.fliplr(np.flipud(matrix))

# Step 5: Map each occurrence of 8 to 9
matrix[matrix == 8] = 9

# Step 6: Remove every 2-nd row (1-indexed)
matrix = np.delete(matrix, 1, axis=0)

# Step 7: Mirror along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix))

# Step 8: Mirror along the diagonal
matrix = matrix.T

# Print the final matrix
print(matrix)