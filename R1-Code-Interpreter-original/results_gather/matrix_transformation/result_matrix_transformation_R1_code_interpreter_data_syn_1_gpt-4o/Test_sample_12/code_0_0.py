def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed = list(zip(*matrix))
    # Reverse each row of the transposed matrix
    rotated = [list(row)[::-1] for row in transposed]
    return rotated

# Original matrix
matrix = [
    ['L', 'U'],
    ['X', 'K']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output as a list of values and symbols
formatted_output = ','.join([','.join(row) for row in rotated_matrix])
print(f"<<<{formatted_output}>>>")