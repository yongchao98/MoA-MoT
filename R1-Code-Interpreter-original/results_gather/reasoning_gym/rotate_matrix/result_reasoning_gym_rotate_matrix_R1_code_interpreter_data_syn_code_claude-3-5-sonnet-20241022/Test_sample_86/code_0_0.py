def rotate_matrix_90(matrix):
    n = len(matrix)
    # Create a new matrix with rotated elements
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    return rotated

# Input matrix
matrix = [
    [4, 1, 5, 4],
    [0, 4, 5, 0],
    [9, 7, 8, 9],
    [5, 9, 8, 2]
]

# Rotate 810 degrees (equivalent to 90 degrees)
result = rotate_matrix_90(matrix)

# Format output
output = "\n".join(" ".join(str(x) for x in row) for row in result)
print(output)