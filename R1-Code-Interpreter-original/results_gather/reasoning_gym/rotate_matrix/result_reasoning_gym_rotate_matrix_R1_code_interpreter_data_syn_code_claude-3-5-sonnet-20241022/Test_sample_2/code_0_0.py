def rotate_90_clockwise(matrix):
    # Transpose
    n = len(matrix)
    for i in range(n):
        for j in range(i, n):
            matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]
    
    # Reverse each row
    for i in range(n):
        matrix[i] = matrix[i][::-1]
    
    return matrix

# Input matrix
matrix = [
    [5, 9, 0, 3, 0, 6, 5, 1, 6, 4],
    [5, 1, 9, 8, 2, 1, 4, 5, 9, 0],
    [2, 8, 1, 0, 5, 4, 8, 0, 6, 2],
    [9, 4, 4, 6, 5, 3, 3, 7, 8, 2],
    [9, 2, 9, 6, 8, 6, 3, 3, 1, 9],
    [1, 7, 4, 6, 1, 9, 4, 3, 2, 6],
    [5, 7, 1, 8, 7, 6, 7, 5, 2, 7],
    [4, 2, 3, 2, 4, 4, 1, 4, 8, 1],
    [9, 8, 4, 6, 3, 6, 4, 8, 5, 9],
    [7, 1, 1, 0, 5, 6, 5, 2, 2, 6]
]

# 450 degrees = 360 + 90 degrees = 1 rotation of 90 degrees
matrix = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in matrix:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())