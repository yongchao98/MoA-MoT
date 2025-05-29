import numpy as np

# Initial matrix
matrix = np.array([
    [0, 5, 7, 8, 1, 7, 9, 4, 6, 6],
    [7, 2, 0, 1, 6, 1, 8, 4, 5, 6],
    [9, 9, 5, 9, 0, 8, 9, 3, 8, 6]
])

# Identity transformation (no change)
identity_matrix = matrix.copy()

# Rotate 270 degrees (90 degrees counterclockwise)
rotated_270 = np.rot90(identity_matrix, k=1)

# Rotate 180 degrees
rotated_180 = np.rot90(rotated_270, k=2)

# Convert the final matrix to the required format
final_matrix = rotated_180.tolist()
formatted_output = '\n'.join(' '.join(map(str, row)) for row in final_matrix)

print(formatted_output)