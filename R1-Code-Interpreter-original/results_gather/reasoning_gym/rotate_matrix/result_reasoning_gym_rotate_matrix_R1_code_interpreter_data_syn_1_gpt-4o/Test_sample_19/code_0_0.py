def rotate_180(matrix):
    # Reverse the order of rows and then reverse each row
    return [row[::-1] for row in matrix[::-1]]

# Given matrix
matrix = [
    [8, 8, 2, 3, 3],
    [8, 4, 6, 8, 7],
    [8, 8, 5, 9, 5],
    [6, 4, 6, 0, 1],
    [0, 8, 4, 0, 7]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)