def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the rotated matrix
    rotated_matrix = [list(reversed(row)) for row in transposed_matrix]
    return rotated_matrix

# Original matrix
matrix = [
    ['N', 'U', 'L'],
    ['M', 'C', 'A'],
    ['Q', 'S', 'T']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output
formatted_output = ','.join([','.join(row) for row in rotated_matrix])
print(f"<<<{formatted_output}>>>")