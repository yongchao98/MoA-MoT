def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the 90-degree rotated matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Input matrix
matrix = [
    [8, 0, 2, 0, 6, 4],
    [5, 0, 5, 2, 6, 4],
    [6, 0, 8, 7, 9, 6],
    [0, 3, 3, 1, 6, 5],
    [3, 2, 6, 3, 5, 8],
    [4, 6, 2, 1, 7, 8]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)