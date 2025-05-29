def rotate_matrix_90_clockwise(matrix):
    n = len(matrix)
    m = len(matrix[0])
    rotated_matrix = [[0] * n for _ in range(m)]
    
    for i in range(n):
        for j in range(m):
            rotated_matrix[j][n - 1 - i] = matrix[i][j]
    
    return rotated_matrix

matrix = [
    [1, 4, 3],
    [6, 8, 7]
]

rotated = rotate_matrix_90_clockwise(matrix)
print(rotated)