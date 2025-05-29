# Input matrix
matrix = [
    [4, 0, 1, 2, 7, 7, 3, 4, 6, 9],
    [6, 9, 1, 7, 1, 5, 8, 5, 7, 0],
    [8, 8, 1, 6, 7, 2, 8, 5, 2, 7],
    [8, 1, 8, 7, 9, 3, 9, 6, 8, 4],
    [5, 6, 0, 6, 0, 6, 3, 5, 9, 1],
    [0, 6, 5, 4, 3, 5, 7, 6, 2, 4],
    [0, 6, 6, 1, 0, 4, 0, 2, 5, 4],
    [0, 7, 0, 2, 4, 8, 7, 0, 9, 9],
    [2, 7, 3, 8, 7, 6, 8, 7, 6, 8],
    [4, 5, 6, 9, 0, 1, 3, 4, 5, 8]
]

# Function to rotate matrix 90 degrees clockwise
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

# Rotate the matrix
rotated = rotate_90_clockwise(matrix)

# Format output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())