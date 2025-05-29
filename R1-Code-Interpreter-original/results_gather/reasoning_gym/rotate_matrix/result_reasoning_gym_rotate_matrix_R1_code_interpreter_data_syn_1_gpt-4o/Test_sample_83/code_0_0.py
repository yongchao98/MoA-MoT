def rotate_270_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse the order of columns
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Input matrix
matrix = [
    [3, 1, 8, 1, 5, 8],
    [4, 4, 5, 7, 0, 8],
    [8, 2, 2, 0, 4, 2],
    [7, 9, 3, 2, 2, 2],
    [8, 3, 2, 6, 0, 3],
    [8, 5, 9, 3, 4, 9]
]

# Rotate the matrix by 270 degrees clockwise
rotated_matrix = rotate_270_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)