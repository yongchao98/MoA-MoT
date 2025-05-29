def rotate_matrix_90_clockwise(matrix):
    n = len(matrix)
    m = len(matrix[0])
    rotated_matrix = [['' for _ in range(n)] for _ in range(m)]
    
    for i in range(n):
        for j in range(m):
            rotated_matrix[j][n-1-i] = matrix[i][j]
    
    return rotated_matrix

matrix = [
    ['1', '8', '4'],
    ['7', '4', '1']
]

rotated = rotate_matrix_90_clockwise(matrix)
result = ','.join([','.join(row) for row in rotated])
print(f"<<<{result}>>>")