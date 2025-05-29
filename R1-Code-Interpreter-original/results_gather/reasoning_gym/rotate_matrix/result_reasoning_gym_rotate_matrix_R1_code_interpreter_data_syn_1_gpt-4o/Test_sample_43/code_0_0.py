def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = [list(row) for row in zip(*matrix)]
    # Reverse each row
    rotated_matrix = [row[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [9, 6, 7, 0, 8, 7, 5, 1, 5],
    [1, 4, 4, 7, 5, 9, 7, 4, 6],
    [0, 6, 8, 9, 7, 3, 2, 3, 2],
    [4, 7, 9, 2, 7, 0, 4, 9, 0],
    [1, 0, 4, 5, 3, 0, 1, 6, 1],
    [4, 5, 8, 1, 7, 4, 4, 9, 3],
    [3, 6, 0, 0, 4, 5, 7, 4, 0],
    [8, 2, 8, 8, 7, 0, 4, 8, 5],
    [2, 5, 4, 6, 7, 7, 1, 6, 4]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)