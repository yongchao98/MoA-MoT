import numpy as np

# Initial matrix
matrix = [
    [8, 7, 8, 5, 7, 0, 0, 2, 5],
    [2, 4, 0, 2, 0, 5, 3, 7, 3],
    [3, 2, 2, 5, 8, 5, 8, 6, 0],
    [5, 4, 2, 1, 8, 6, 6, 7, 3],
    [2, 0, 7, 7, 3, 5, 9, 3, 2]
]

# Step 1: Identity transformation (no change)
# matrix remains the same

# Step 2: Rotate the matrix 90 degrees clockwise
rotated_matrix = np.array(matrix).T[:, ::-1]

# Step 3: Remove every 2-nd row (1-indexed)
final_matrix = np.delete(rotated_matrix, np.arange(1, rotated_matrix.shape[0], 2), axis=0)

# Convert the final matrix to a list of lists and print it
final_matrix_list = final_matrix.tolist()
for row in final_matrix_list:
    print(' '.join(map(str, row)))