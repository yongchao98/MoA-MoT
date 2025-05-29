def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the 90-degree rotated matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Given matrix
matrix = [
    [9, 1, 5, 3, 3],
    [7, 5, 9, 3, 6],
    [8, 8, 3, 9, 9],
    [9, 7, 9, 6, 1],
    [0, 8, 1, 6, 7]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)