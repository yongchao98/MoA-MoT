def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed = list(zip(*matrix))
    # Reverse each row
    rotated = [list(row)[::-1] for row in transposed]
    return rotated

# Initial matrix
matrix = [
    [4, 2],
    [7, 5]
]

# Rotate 90 degrees clockwise twice
for _ in range(2):
    matrix = rotate_90_clockwise(matrix)

print(matrix)