def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated_matrix = [list(row)[::-1] for row in transposed_matrix]
    return rotated_matrix

# Define the matrix
matrix = [
    [4, 2, 7],
    [7, 8, 6]
]

# Rotate the matrix
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the result
print(rotated_matrix)