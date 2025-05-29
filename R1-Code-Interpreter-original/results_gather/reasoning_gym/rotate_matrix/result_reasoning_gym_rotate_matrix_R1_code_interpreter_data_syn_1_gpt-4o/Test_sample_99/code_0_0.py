def rotate_270_clockwise(matrix):
    # Transpose the matrix
    transposed = list(zip(*matrix))
    # Reverse the order of columns
    rotated = [list(row)[::-1] for row in transposed]
    return rotated

# Given matrix
matrix = [
    [5, 8, 1, 6, 4, 0, 6],
    [1, 5, 9, 9, 1, 2, 3],
    [8, 7, 7, 3, 3, 8, 0],
    [2, 7, 9, 3, 6, 1, 0],
    [7, 4, 6, 1, 2, 1, 2],
    [1, 9, 5, 1, 8, 9, 6],
    [3, 2, 9, 7, 5, 0, 2]
]

# Rotate the matrix
rotated_matrix = rotate_270_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)