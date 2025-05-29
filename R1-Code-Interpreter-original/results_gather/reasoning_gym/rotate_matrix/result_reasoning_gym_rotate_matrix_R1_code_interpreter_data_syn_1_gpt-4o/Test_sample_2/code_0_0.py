def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Input matrix
matrix = [
    [5, 9, 0, 3, 0, 6, 5, 1, 6, 4],
    [5, 1, 9, 8, 2, 1, 4, 5, 9, 0],
    [2, 8, 1, 0, 5, 4, 8, 0, 6, 2],
    [9, 4, 4, 6, 5, 3, 3, 7, 8, 2],
    [9, 2, 9, 6, 8, 6, 3, 3, 1, 9],
    [1, 7, 4, 6, 1, 9, 4, 3, 2, 6],
    [5, 7, 1, 8, 7, 6, 7, 5, 2, 7],
    [4, 2, 3, 2, 4, 4, 1, 4, 8, 1],
    [9, 8, 4, 6, 3, 6, 4, 8, 5, 9],
    [7, 1, 1, 0, 5, 6, 5, 2, 2, 6]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)