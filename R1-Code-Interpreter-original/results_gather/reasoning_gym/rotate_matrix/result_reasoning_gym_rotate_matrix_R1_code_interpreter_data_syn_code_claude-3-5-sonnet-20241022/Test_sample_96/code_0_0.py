# Input matrix
matrix = [
    [6, 8, 5, 2, 7],
    [0, 2, 5, 8, 4],
    [9, 8, 4, 1, 1],
    [0, 4, 5, 1, 7],
    [2, 6, 8, 5, 7]
]

def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Transpose matrix
    for i in range(n):
        for j in range(i, n):
            matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]
    
    # Reverse each row
    for i in range(n):
        matrix[i] = matrix[i][::-1]
    
    return matrix

# Rotate matrix
rotated = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())