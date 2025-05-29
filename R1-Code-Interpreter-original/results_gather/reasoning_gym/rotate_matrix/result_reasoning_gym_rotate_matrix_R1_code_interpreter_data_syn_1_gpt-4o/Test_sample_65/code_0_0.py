def rotate_matrix_270_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse the order of the rows
    rotated_matrix = transposed_matrix[::-1]
    # Convert tuples back to lists
    rotated_matrix = [list(row) for row in rotated_matrix]
    return rotated_matrix

# Original matrix
matrix = [
    [6, 5, 0, 9, 0],
    [8, 0, 0, 4, 5],
    [0, 1, 7, 7, 1],
    [5, 3, 4, 2, 4],
    [6, 8, 3, 2, 1]
]

# Rotate the matrix by 270 degrees clockwise
rotated_matrix = rotate_matrix_270_clockwise(matrix)
print(rotated_matrix)