def rotate_90_clockwise(matrix):
    # Transpose the matrix and then reverse each row
    return [list(reversed(col)) for col in zip(*matrix)]

# Original matrix
matrix = [
    ['G', 'Y', 'Z'],
    ['A', 'A', 'A'],
    ['U', 'B', 'Y']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output as a single list of values
output = [item for row in rotated_matrix for item in row]
print(output)