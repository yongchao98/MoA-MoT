import numpy as np

# Initial matrix
matrix = np.array([
    [8, 6, 5, 5, 8, 6, 5, 0],
    [8, 6, 3, 9, 0, 6, 3, 2],
    [8, 9, 0, 2, 7, 9, 2, 2],
    [3, 5, 5, 9, 2, 3, 4, 4],
    [1, 5, 1, 0, 0, 2, 8, 4],
    [4, 5, 9, 9, 0, 5, 7, 1]
])

# Step 2: Set all elements divisible by 9 to zero
matrix[matrix % 9 == 0] = 0

# Step 3: Mirror the matrix along the diagonal (transpose)
matrix = matrix.T

# Step 4: Crop the matrix to rows 8-8 and columns 3-5 (1-indexed)
# Adjusting for 0-indexing, this means row 5 and columns 2-4
cropped_matrix = matrix[2:5, 2:5]

# Step 5: Remove every 3rd column (1-indexed)
# In this case, we remove the 3rd column of the cropped matrix
final_matrix = np.delete(cropped_matrix, 2, axis=1)

# Print the final matrix
print(final_matrix)