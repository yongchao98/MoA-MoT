def rotate_180(matrix):
    # Reverse each row
    matrix = [row[::-1] for row in matrix]
    # Reverse the order of rows
    matrix.reverse()
    return matrix

# Input matrix
matrix = [
    [3, 9, 4, 1, 5, 3, 7, 9, 2, 1],
    [4, 7, 3, 8, 6, 6, 2, 5, 4, 7],
    [1, 6, 1, 8, 7, 9, 1, 8, 9, 7],
    [4, 5, 9, 0, 8, 5, 6, 6, 2, 5],
    [6, 1, 7, 1, 9, 0, 0, 3, 7, 6],
    [4, 5, 7, 9, 3, 0, 3, 9, 5, 7],
    [1, 7, 5, 9, 1, 3, 0, 7, 1, 5],
    [1, 4, 4, 9, 8, 3, 6, 9, 9, 3],
    [3, 2, 2, 6, 5, 5, 7, 2, 6, 5],
    [0, 6, 1, 8, 9, 6, 2, 0, 7, 2]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)