def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the 90-degree clockwise rotation
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [5, 6, 3, 1, 5, 0, 0],
    [2, 6, 5, 1, 4, 5, 5],
    [4, 2, 8, 0, 9, 7, 8],
    [7, 7, 5, 3, 5, 8, 8],
    [5, 0, 2, 0, 1, 7, 8],
    [3, 5, 1, 2, 9, 1, 8],
    [4, 9, 4, 4, 3, 6, 6]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)