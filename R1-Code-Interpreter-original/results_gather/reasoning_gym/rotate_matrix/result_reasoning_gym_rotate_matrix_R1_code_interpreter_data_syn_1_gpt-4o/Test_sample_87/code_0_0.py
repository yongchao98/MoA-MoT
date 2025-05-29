def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the 90-degree clockwise rotated matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [2, 8, 8, 6, 1, 1],
    [6, 3, 2, 4, 4, 3],
    [7, 9, 0, 6, 4, 6],
    [0, 5, 3, 2, 7, 0],
    [8, 8, 7, 6, 1, 0],
    [3, 3, 7, 6, 7, 8]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)