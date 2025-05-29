def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed = list(zip(*matrix))
    # Reverse each row to get the 90 degree rotated matrix
    rotated = [list(reversed(row)) for row in transposed]
    return rotated

# Original matrix
matrix = [
    ['Z', 'P', 'E'],
    ['Z', 'T', 'V'],
    ['G', 'I', 'V']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output
formatted_output = ','.join([','.join(row) for row in rotated_matrix])
print(formatted_output)