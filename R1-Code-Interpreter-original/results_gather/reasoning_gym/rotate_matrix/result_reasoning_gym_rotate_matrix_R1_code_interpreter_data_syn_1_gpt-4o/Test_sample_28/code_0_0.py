def rotate_270_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = [list(row) for row in zip(*matrix)]
    # Reverse the order of columns
    rotated_matrix = [row[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [4, 9, 2, 4, 3, 7],
    [7, 5, 0, 1, 2, 2],
    [6, 0, 3, 0, 6, 2],
    [2, 4, 9, 2, 8, 9],
    [8, 9, 0, 9, 9, 1],
    [4, 2, 0, 9, 3, 4]
]

# Rotate the matrix by 270 degrees clockwise
rotated_matrix = rotate_270_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)