def rotate_180(matrix):
    # Reverse the order of rows and then reverse each row
    return [row[::-1] for row in matrix[::-1]]

# Given matrix
matrix = [
    [5, 3, 2, 6, 9, 4, 8, 4],
    [6, 9, 0, 3, 0, 7, 2, 2],
    [5, 2, 4, 9, 9, 0, 1, 9],
    [0, 6, 4, 3, 1, 5, 7, 2],
    [8, 0, 1, 1, 2, 2, 0, 4],
    [5, 3, 9, 0, 9, 9, 2, 6],
    [0, 4, 6, 2, 1, 1, 6, 0],
    [7, 6, 3, 0, 1, 2, 1, 7]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)