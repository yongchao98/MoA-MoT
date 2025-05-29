def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated_matrix = [list(reversed(row)) for row in transposed_matrix]
    return rotated_matrix

# Original matrix
matrix = [
    ['7', '1', '2'],
    ['7', '8', '9']
]

# Rotate the matrix
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Format the output as a list of strings
formatted_output = [','.join(row) for row in rotated_matrix]
print(formatted_output)