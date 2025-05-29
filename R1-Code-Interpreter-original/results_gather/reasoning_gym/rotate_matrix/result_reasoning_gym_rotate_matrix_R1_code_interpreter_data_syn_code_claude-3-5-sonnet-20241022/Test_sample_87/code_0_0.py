def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Transpose the matrix
    for i in range(n):
        for j in range(i, n):
            matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]
    
    # Reverse each row
    for i in range(n):
        matrix[i] = matrix[i][::-1]
    return matrix

# Input matrix
matrix = [
    [2, 8, 8, 6, 1, 1],
    [6, 3, 2, 4, 4, 3],
    [7, 9, 0, 6, 4, 6],
    [0, 5, 3, 2, 7, 0],
    [8, 8, 7, 6, 1, 0],
    [3, 3, 7, 6, 7, 8]
]

# Rotate once (90 degrees clockwise)
rotated = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())