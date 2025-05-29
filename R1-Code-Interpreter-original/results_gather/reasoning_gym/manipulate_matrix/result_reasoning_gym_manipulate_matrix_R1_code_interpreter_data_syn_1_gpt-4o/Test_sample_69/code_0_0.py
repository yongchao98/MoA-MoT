import numpy as np

# Initial matrix
matrix = np.array([
    [5, 7, 0, 8, 3, 9, 3, 4, 7, 3],
    [5, 7, 1, 4, 3, 2, 1, 7, 8, 5],
    [3, 7, 5, 3, 6, 6, 3, 1, 6, 7],
    [4, 3, 1, 9, 5, 4, 6, 6, 3, 6],
    [1, 6, 6, 0, 7, 3, 8, 1, 5, 0]
])

# Step 1: Identity transformation (no change)
# Step 2: Crop to rows 1-3 and columns 10-10 (1-indexed)
cropped_matrix = matrix[0:3, 9:10]

# Step 3: Mirror along the counterdiagonal
mirrored_counterdiagonal = np.fliplr(np.flipud(cropped_matrix.T))

# Step 4: Vertically mirror the matrix
vertically_mirrored = np.flipud(mirrored_counterdiagonal)

# Step 5: Crop to rows 1-1 and columns 3-3 (1-indexed)
# Since the matrix is 1 column, this will result in an empty operation
cropped_matrix_2 = vertically_mirrored[0:1, 2:3]

# Step 6: Horizontally mirror the matrix
horizontally_mirrored = np.fliplr(cropped_matrix_2)

# Step 7: Map each occurrence of 2 to 8
final_matrix = np.where(horizontally_mirrored == 2, 8, horizontally_mirrored)

# Output the final matrix
print(final_matrix)