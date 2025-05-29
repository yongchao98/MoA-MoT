def rotate_90_clockwise(matrix):
    return [list(reversed(col)) for col in zip(*matrix)]

# Initial matrix
matrix = [
    [8, 5],
    [0, 1]
]

# Rotate 90 degrees clockwise four times (equivalent to 720 degrees)
for _ in range(4):
    matrix = rotate_90_clockwise(matrix)

print(matrix)