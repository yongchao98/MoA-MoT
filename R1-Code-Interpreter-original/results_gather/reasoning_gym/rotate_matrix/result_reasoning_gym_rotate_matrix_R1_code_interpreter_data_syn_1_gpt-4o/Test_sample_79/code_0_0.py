def rotate_180(matrix):
    # Reverse each row
    matrix = [row[::-1] for row in matrix]
    # Reverse the order of rows
    matrix.reverse()
    return matrix

# Given matrix
matrix = [
    [5, 7, 3, 6, 6],
    [9, 2, 3, 6, 8],
    [3, 2, 7, 6, 6],
    [4, 6, 3, 5, 7],
    [0, 5, 3, 9, 7]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)