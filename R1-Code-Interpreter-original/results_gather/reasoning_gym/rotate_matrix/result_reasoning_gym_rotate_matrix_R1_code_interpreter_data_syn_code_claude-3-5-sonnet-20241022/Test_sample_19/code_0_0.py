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
    [8, 8, 2, 3, 3],
    [8, 4, 6, 8, 7],
    [8, 8, 5, 9, 5],
    [6, 4, 6, 0, 1],
    [0, 8, 4, 0, 7]
]

# For 900 degrees = 2 complete rotations (720°) + half rotation (180°)
# We need to rotate 2 times

# First rotation (90°)
matrix = rotate_90_clockwise(matrix)
# Second rotation (180°)
matrix = rotate_90_clockwise(matrix)
# Third rotation (270°)
matrix = rotate_90_clockwise(matrix)
# Fourth rotation (360°)
matrix = rotate_90_clockwise(matrix)
# Fifth rotation (450°)
matrix = rotate_90_clockwise(matrix)
# Sixth rotation (540°)
matrix = rotate_90_clockwise(matrix)
# Seventh rotation (630°)
matrix = rotate_90_clockwise(matrix)
# Eighth rotation (720°)
matrix = rotate_90_clockwise(matrix)
# Ninth rotation (810°)
matrix = rotate_90_clockwise(matrix)
# Tenth rotation (900°)
matrix = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in matrix:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())