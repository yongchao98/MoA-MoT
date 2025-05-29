def rotate_90_clockwise(matrix):
    # Transpose the matrix and then reverse each row
    return [''.join(row) for row in zip(*matrix[::-1])]

# Original matrix
matrix = [
    "FKI",
    "FQY",
    "CIE"
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the result
print(rotated_matrix)