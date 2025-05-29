import numpy as np

# Initial matrix
matrix = np.array([
    [7, 2, 6, 5, 6, 8, 1, 8, 1],
    [1, 7, 4, 0, 3, 5, 3, 4, 3],
    [2, 2, 3, 5, 5, 6, 3, 6, 7],
    [8, 4, 6, 2, 1, 2, 0, 9, 6]
])

# Step 2: Crop the matrix to rows 2-2 and columns 3-9 (1-indexed)
cropped_matrix = matrix[1:2, 2:9]

# Step 3: Rotate the matrix 180 degrees
rotated_matrix = np.rot90(cropped_matrix, 2)

# Step 4: Remove every 7-th column (1-indexed)
# Since the cropped matrix has only one row, we need to check if it has 7 columns
# In this case, it doesn't, so this step will not remove any column
final_matrix = rotated_matrix

# Print the final matrix
print(final_matrix)