def rotate_matrix_360(matrix):
    # A 360-degree rotation results in the same matrix
    return matrix

# Define the matrix
matrix = [
    [1, 8, 7, 9, 2, 7],
    [4, 9, 1, 7, 6, 8],
    [8, 3, 7, 1, 0, 4],
    [6, 0, 6, 5, 8, 6],
    [7, 8, 6, 2, 0, 4],
    [9, 9, 3, 1, 7, 5]
]

# Rotate the matrix by 360 degrees
rotated_matrix = rotate_matrix_360(matrix)

# Print the rotated matrix
print(rotated_matrix)