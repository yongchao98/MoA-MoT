def rotate_270_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = [list(row) for row in zip(*matrix)]
    # Reverse the order of columns in the transposed matrix
    rotated_matrix = [row[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [3, 3, 9, 7, 1, 5, 3, 3, 4, 5],
    [5, 2, 1, 5, 4, 2, 4, 4, 4, 3],
    [3, 8, 8, 1, 6, 3, 3, 0, 1, 0],
    [7, 6, 6, 9, 6, 6, 8, 6, 1, 0],
    [0, 6, 5, 4, 9, 3, 2, 6, 6, 9],
    [5, 7, 4, 5, 6, 4, 4, 6, 4, 8],
    [8, 5, 0, 5, 8, 3, 8, 3, 5, 5],
    [8, 4, 2, 8, 9, 1, 0, 5, 3, 9],
    [1, 3, 2, 6, 0, 9, 2, 7, 1, 7],
    [0, 1, 9, 3, 2, 8, 5, 5, 3, 7]
]

# Rotate the matrix
rotated_matrix = rotate_270_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)