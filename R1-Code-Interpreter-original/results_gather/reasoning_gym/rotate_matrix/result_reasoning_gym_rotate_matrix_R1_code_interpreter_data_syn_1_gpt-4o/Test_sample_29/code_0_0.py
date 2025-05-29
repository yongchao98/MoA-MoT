def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [2, 8, 3, 5, 4, 5, 0, 7, 5],
    [5, 1, 6, 7, 5, 2, 5, 0, 8],
    [1, 2, 5, 0, 6, 2, 2, 6, 4],
    [4, 5, 7, 5, 9, 5, 2, 9, 4],
    [6, 6, 3, 2, 2, 0, 8, 5, 7],
    [0, 8, 6, 1, 4, 0, 6, 2, 1],
    [7, 5, 1, 6, 3, 1, 9, 5, 9],
    [5, 0, 7, 4, 6, 6, 3, 2, 1],
    [6, 7, 1, 4, 7, 6, 0, 5, 3]
]

# Rotate the matrix 90 degrees clockwise
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)