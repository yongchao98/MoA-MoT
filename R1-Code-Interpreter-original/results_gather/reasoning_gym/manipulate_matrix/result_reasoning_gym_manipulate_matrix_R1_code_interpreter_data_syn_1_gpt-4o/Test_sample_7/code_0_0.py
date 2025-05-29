import numpy as np

# Initial matrix
matrix = np.array([
    [6, 0, 3, 3, 6, 9],
    [0, 4, 5, 2, 0, 4],
    [1, 5, 0, 1, 6, 5],
    [6, 1, 6, 8, 9, 4],
    [4, 1, 7, 6, 4, 9]
])

# Step 2: Map each occurrence of 0 to 3
matrix[matrix == 0] = 3

# Step 3: Remove every 6-th column (1-indexed)
matrix = np.delete(matrix, 5, axis=1)

# Step 4: Mirror the matrix along the counterdiagonal
matrix = np.fliplr(np.flipud(matrix.T))

# Step 5: Crop the matrix to rows 4-5 and columns 3-4 (1-indexed)
cropped_matrix = matrix[3:5, 2:4]

# Step 6: Set all elements divisible by 8 to zero
cropped_matrix[cropped_matrix % 8 == 0] = 0

# Step 7: Horizontally mirror the matrix
final_matrix = np.fliplr(cropped_matrix)

# Print the final matrix
print(final_matrix)