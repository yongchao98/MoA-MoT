def rotate_matrix_270_clockwise(matrix):
    # Transpose the matrix and then reverse each row to achieve 270-degree clockwise rotation
    return [list(row) for row in zip(*matrix[::-1])]

# Input matrix
matrix = [
    [7, 6, 3, 4, 0, 6, 3],
    [4, 4, 7, 7, 3, 1, 5],
    [4, 3, 5, 8, 0, 0, 6],
    [9, 9, 6, 3, 5, 8, 7],
    [0, 1, 6, 2, 1, 3, 1],
    [4, 9, 1, 8, 1, 0, 0],
    [3, 7, 7, 2, 9, 9, 6]
]

# Rotate the matrix
rotated_matrix = rotate_matrix_270_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)