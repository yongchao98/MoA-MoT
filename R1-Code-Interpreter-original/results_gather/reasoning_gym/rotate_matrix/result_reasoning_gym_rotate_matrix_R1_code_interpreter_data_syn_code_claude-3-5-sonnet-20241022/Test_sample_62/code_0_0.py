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
    [2, 1, 4, 4, 6, 0, 7, 5, 4, 6],
    [9, 4, 8, 0, 6, 0, 6, 5, 7, 5],
    [3, 0, 6, 7, 2, 4, 2, 1, 4, 7],
    [4, 1, 5, 6, 5, 0, 8, 0, 4, 3],
    [0, 7, 5, 3, 4, 2, 0, 2, 7, 9],
    [3, 7, 0, 3, 8, 8, 1, 9, 9, 4],
    [0, 2, 8, 2, 3, 3, 0, 6, 9, 9],
    [9, 2, 3, 2, 3, 7, 2, 8, 6, 3],
    [5, 5, 6, 2, 8, 9, 2, 9, 5, 6],
    [2, 7, 2, 8, 5, 5, 4, 8, 0, 5]
]

# 450 degrees = 360 + 90 degrees = 1 rotation of 90 degrees
result = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())