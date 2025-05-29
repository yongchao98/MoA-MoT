def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Transpose
    for i in range(n):
        for j in range(i, n):
            matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]
    
    # Reverse each row
    for i in range(n):
        matrix[i] = matrix[i][::-1]
    return matrix

# Input matrix
matrix = [
    [9, 6, 7, 0, 8, 7, 5, 1, 5],
    [1, 4, 4, 7, 5, 9, 7, 4, 6],
    [0, 6, 8, 9, 7, 3, 2, 3, 2],
    [4, 7, 9, 2, 7, 0, 4, 9, 0],
    [1, 0, 4, 5, 3, 0, 1, 6, 1],
    [4, 5, 8, 1, 7, 4, 4, 9, 3],
    [3, 6, 0, 0, 4, 5, 7, 4, 0],
    [8, 2, 8, 8, 7, 0, 4, 8, 5],
    [2, 5, 4, 6, 7, 7, 1, 6, 4]
]

# Rotate once (90 degrees)
result = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())