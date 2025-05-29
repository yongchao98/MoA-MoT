import numpy as np

# Initial matrix
matrix = np.array([
    [7, 7, 8, 0, 1],
    [3, 8, 7, 4, 9],
    [8, 1, 3, 7, 8],
    [4, 6, 2, 4, 8],
    [8, 5, 1, 8, 2],
    [0, 5, 9, 3, 2]
])

# Step 2: Crop the matrix to rows 4-6 and columns 4-4 (1-indexed)
cropped_matrix = matrix[3:6, 3:4]

# Step 3: Mirror the matrix along the counterdiagonal
mirrored_counterdiagonal = np.fliplr(np.flipud(cropped_matrix.T))

# Step 4: Horizontally mirror the matrix
horizontally_mirrored = np.fliplr(mirrored_counterdiagonal)

# Step 5: Rotate the matrix 90 degrees
rotated_90 = np.rot90(horizontally_mirrored, -1)

# Step 6: Mirror the matrix along the diagonal (transpose)
mirrored_diagonal = rotated_90.T

# Step 7: Mirror the matrix along the diagonal (transpose again)
mirrored_diagonal_again = mirrored_diagonal.T

# Step 8: Mirror the matrix along the counterdiagonal
mirrored_counterdiagonal_again = np.fliplr(np.flipud(mirrored_diagonal_again.T))

# Step 9: Remove every 3rd column (1-indexed)
# Since the matrix is 1 column wide, this step has no effect

# Step 10: Set all elements divisible by 8 to zero
mirrored_counterdiagonal_again[mirrored_counterdiagonal_again % 8 == 0] = 0

# Step 11: Map each occurrence of 1 to 0
mirrored_counterdiagonal_again[mirrored_counterdiagonal_again == 1] = 0

# Print the final matrix
print(mirrored_counterdiagonal_again)