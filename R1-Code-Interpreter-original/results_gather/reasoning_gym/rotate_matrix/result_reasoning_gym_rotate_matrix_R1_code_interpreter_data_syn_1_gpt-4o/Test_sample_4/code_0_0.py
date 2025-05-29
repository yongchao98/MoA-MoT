def rotate_270_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the 270-degree clockwise rotation
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Define the input matrix
matrix = [
    [0, 0, 4, 1, 0, 8, 2],
    [3, 4, 6, 5, 8, 4, 8],
    [1, 9, 6, 7, 6, 5, 6],
    [3, 6, 2, 1, 8, 7, 3],
    [8, 4, 2, 7, 9, 5, 8],
    [2, 0, 9, 6, 7, 2, 3],
    [6, 5, 3, 7, 4, 9, 0]
]

# Rotate the matrix
rotated_matrix = rotate_270_clockwise(matrix)

# Print the rotated matrix
for row in rotated_matrix:
    print(' '.join(map(str, row)))