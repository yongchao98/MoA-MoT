import numpy as np

# Initial matrix
matrix = np.array([
    [3, 9, 7, 1, 5],
    [3, 3, 4, 5, 5],
    [2, 1, 5, 4, 2],
    [4, 4, 4, 3, 3],
    [8, 8, 1, 6, 3],
    [3, 0, 1, 0, 7],
    [6, 6, 9, 6, 6],
    [8, 6, 1, 0, 0],
    [6, 5, 4, 9, 3],
    [2, 6, 6, 9, 5]
])

# Step 2: Crop the matrix to rows 6-9 and columns 3-4 (1-indexed)
cropped_matrix = matrix[5:9, 2:4]

# Step 3: Remove every 2-nd column (1-indexed)
cropped_matrix = np.delete(cropped_matrix, 1, axis=1)

# Step 4: Vertically mirror the matrix
mirrored_matrix = np.flipud(cropped_matrix)

# Step 5: Remove every 4-th row (1-indexed)
mirrored_matrix = np.delete(mirrored_matrix, 3, axis=0)

# Step 6: Mirror the matrix along the diagonal (transpose)
transposed_matrix = mirrored_matrix.T

# Step 7: Remove every 2-nd column (1-indexed)
transposed_matrix = np.delete(transposed_matrix, 1, axis=1)

# Step 8: Map each occurrence of 3 to 5
transposed_matrix[transposed_matrix == 3] = 5

# Step 9: Mirror the matrix along the diagonal (transpose)
final_matrix = transposed_matrix.T

# Print the final matrix
print(final_matrix)