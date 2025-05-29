def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed = list(zip(*matrix))
    # Reverse each row to get the rotated matrix
    rotated = [list(reversed(row)) for row in transposed]
    return rotated

# Original matrix
matrix = [
    ['C', 'T', 'K'],
    ['Y', 'H', 'F'],
    ['O', 'N', 'O']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output as a single list of values
output = []
for row in rotated_matrix:
    output.extend(row)

print(output)