def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Original matrix
matrix = [
    [6, 8, 5, 2, 7],
    [0, 2, 5, 8, 4],
    [9, 8, 4, 1, 1],
    [0, 4, 5, 1, 7],
    [2, 6, 8, 5, 7]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)