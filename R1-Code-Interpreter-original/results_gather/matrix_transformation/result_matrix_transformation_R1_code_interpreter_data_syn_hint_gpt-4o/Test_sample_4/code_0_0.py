def rotate_90_clockwise(matrix):
    # Transpose the matrix and then reverse each row
    return [list(reversed(col)) for col in zip(*matrix)]

# Original matrix
matrix = [
    ['Z', 'T'],
    ['O', 'H']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)