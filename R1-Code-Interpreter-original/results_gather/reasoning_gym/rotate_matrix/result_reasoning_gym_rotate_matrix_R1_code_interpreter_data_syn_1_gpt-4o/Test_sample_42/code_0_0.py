def rotate_180(matrix):
    # Reverse the order of rows and then reverse each row
    return [row[::-1] for row in matrix[::-1]]

# Input matrix
matrix = [
    [1, 5, 1, 4, 2, 3, 5, 1, 7, 7],
    [8, 4, 7, 2, 5, 4, 8, 9, 0, 8],
    [0, 2, 3, 7, 2, 9, 6, 3, 4, 4],
    [6, 0, 8, 3, 2, 9, 0, 3, 5, 0],
    [0, 0, 1, 6, 9, 7, 6, 2, 2, 9],
    [2, 6, 7, 7, 6, 9, 1, 3, 2, 0],
    [7, 7, 2, 7, 2, 3, 4, 0, 4, 8],
    [4, 8, 4, 5, 2, 5, 1, 3, 9, 2],
    [8, 5, 6, 3, 5, 1, 1, 9, 7, 3],
    [2, 9, 4, 4, 8, 2, 0, 8, 9, 7]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)