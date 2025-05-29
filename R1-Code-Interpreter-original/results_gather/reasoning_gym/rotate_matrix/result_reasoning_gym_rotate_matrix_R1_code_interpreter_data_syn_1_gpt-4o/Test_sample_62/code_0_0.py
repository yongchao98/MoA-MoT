def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Input matrix
matrix = [
    [2, 1, 4, 4, 6, 0, 7, 5, 4, 6],
    [9, 4, 8, 0, 6, 0, 6, 5, 7, 5],
    [3, 0, 6, 7, 2, 4, 2, 1, 4, 7],
    [4, 1, 5, 6, 5, 0, 8, 0, 4, 3],
    [0, 7, 5, 3, 4, 2, 0, 2, 7, 9],
    [3, 7, 0, 3, 8, 8, 1, 9, 9, 4],
    [0, 2, 8, 2, 3, 3, 0, 6, 9, 9],
    [9, 2, 3, 2, 3, 7, 2, 8, 6, 3],
    [5, 5, 6, 2, 8, 9, 2, 9, 5, 6],
    [2, 7, 2, 8, 5, 5, 4, 8, 0, 5]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
for row in rotated_matrix:
    print(' '.join(map(str, row)))