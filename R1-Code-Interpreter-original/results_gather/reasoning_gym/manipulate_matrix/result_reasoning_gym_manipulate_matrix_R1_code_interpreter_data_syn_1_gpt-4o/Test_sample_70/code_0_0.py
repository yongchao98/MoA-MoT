import numpy as np

# Original matrix
matrix = np.array([
    [2, 8, 9, 2, 6, 4, 9],
    [6, 9, 9, 0, 8, 8, 9],
    [8, 2, 6, 7, 3, 8, 2],
    [1, 0, 2, 9, 9, 7, 8]
])

# Step 2: Crop the matrix to rows 2-3 and columns 5-6 (1-indexed)
cropped_matrix = matrix[1:3, 4:6]

# Step 3: Set all elements divisible by 2 to zero
cropped_matrix[cropped_matrix % 2 == 0] = 0

# Step 4: Vertically mirror the matrix
mirrored_matrix = np.flipud(cropped_matrix)

# Step 5: Rotate the matrix 270 degrees (90 degrees counterclockwise)
final_matrix = np.rot90(mirrored_matrix, k=1)

# Print the final matrix
print(final_matrix)