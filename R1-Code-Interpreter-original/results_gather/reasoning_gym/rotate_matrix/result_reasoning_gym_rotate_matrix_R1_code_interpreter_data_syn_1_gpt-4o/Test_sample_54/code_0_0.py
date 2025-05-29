def rotate_90_clockwise(matrix):
    return [list(reversed(col)) for col in zip(*matrix)]

# Original matrix
matrix = [
    [8, 5],
    [0, 1]
]

# Rotate 90 degrees clockwise twice (equivalent to 720 degrees)
rotated_matrix = rotate_90_clockwise(matrix)
rotated_matrix = rotate_90_clockwise(rotated_matrix)

print(rotated_matrix)