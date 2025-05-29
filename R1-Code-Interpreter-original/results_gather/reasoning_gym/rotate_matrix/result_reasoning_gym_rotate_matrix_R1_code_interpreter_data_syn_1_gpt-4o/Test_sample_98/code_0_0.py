def rotate_270_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse the order of columns
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Input matrix
matrix = [
    [6, 4, 7, 7, 0, 5, 3, 9, 5, 3],
    [5, 8, 1, 8, 3, 5, 0, 6, 0, 5],
    [2, 0, 7, 6, 5, 7, 6, 6, 1, 0],
    [8, 8, 8, 8, 2, 8, 2, 9, 4, 1],
    [2, 2, 1, 7, 6, 1, 4, 6, 6, 1],
    [5, 8, 2, 2, 5, 1, 4, 6, 8, 9],
    [6, 2, 6, 4, 6, 3, 3, 2, 2, 8],
    [5, 8, 7, 9, 3, 8, 2, 9, 4, 0],
    [5, 8, 0, 0, 9, 6, 5, 7, 7, 8],
    [1, 9, 1, 4, 0, 5, 7, 4, 2, 5]
]

# Rotate the matrix
rotated_matrix = rotate_270_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)