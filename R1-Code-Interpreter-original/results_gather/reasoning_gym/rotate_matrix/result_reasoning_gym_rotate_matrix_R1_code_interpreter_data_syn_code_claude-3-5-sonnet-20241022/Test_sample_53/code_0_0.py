def rotate_matrix_90(matrix):
    n = len(matrix)
    # Create a new matrix with rotated values
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    return rotated

# Input matrix
matrix = [
    [9, 1, 5, 3, 3],
    [7, 5, 9, 3, 6],
    [8, 8, 3, 9, 9],
    [9, 7, 9, 6, 1],
    [0, 8, 1, 6, 7]
]

# Since 810° = 2 complete rotations + 90°
# We only need to rotate once by 90°
result = rotate_matrix_90(matrix)

# Format the output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())