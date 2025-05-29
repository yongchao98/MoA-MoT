def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [4, 1, 5, 4],
    [0, 4, 5, 0],
    [9, 7, 8, 9],
    [5, 9, 8, 2]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)