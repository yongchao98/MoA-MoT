def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Transpose
    for i in range(n):
        for j in range(i, n):
            matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]
    # Reverse rows
    for i in range(n):
        matrix[i] = matrix[i][::-1]
    return matrix

# Input matrix
matrix = [
    [8, 5, 7, 0, 8],
    [3, 9, 3, 4, 7],
    [3, 5, 7, 1, 4],
    [3, 2, 1, 7, 8],
    [5, 3, 7, 5, 3]
]

# 900 degrees = 2 full rotations (720°) + half rotation (180°)
# For 180° rotation, we need two 90° rotations

# Perform rotations
for _ in range(6):  # 2 full rotations (4 times) + half rotation (2 times) = 6 times
    matrix = rotate_90_clockwise(matrix)

# Format output
result = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(result)