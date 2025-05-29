def rotate_90_clockwise(matrix):
    n = len(matrix)
    result = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    return result

# Input matrix
matrix = [
    [7, 6, 3, 4, 0, 6, 3],
    [4, 4, 7, 7, 3, 1, 5],
    [4, 3, 5, 8, 0, 0, 6],
    [9, 9, 6, 3, 5, 8, 7],
    [0, 1, 6, 2, 1, 3, 1],
    [4, 9, 1, 8, 1, 0, 0],
    [3, 7, 7, 2, 9, 9, 6]
]

# Rotate 3 times (270 degrees)
for _ in range(3):
    matrix = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in matrix:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())