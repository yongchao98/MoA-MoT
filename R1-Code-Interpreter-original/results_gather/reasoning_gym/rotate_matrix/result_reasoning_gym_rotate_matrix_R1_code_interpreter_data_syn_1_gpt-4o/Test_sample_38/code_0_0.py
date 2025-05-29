def rotate_180(matrix):
    # Reverse each row
    matrix = [row[::-1] for row in matrix]
    # Reverse the order of rows
    matrix.reverse()
    return matrix

# Original matrix
matrix = [
    [4, 7, 7],
    [2, 5, 1],
    [2, 9, 3]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)