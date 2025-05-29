def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Transpose the matrix
    for i in range(n):
        for j in range(i, n):
            matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]
    
    # Reverse each row
    for i in range(n):
        matrix[i] = matrix[i][::-1]
    
    return matrix

# Input matrix
matrix_str = """6 4 7 7 0 5 3 9 5 3
5 8 1 8 3 5 0 6 0 5
2 0 7 6 5 7 6 6 1 0
8 8 8 8 2 8 2 9 4 1
2 2 1 7 6 1 4 6 6 1
5 8 2 2 5 1 4 6 8 9
6 2 6 4 6 3 3 2 2 8
5 8 7 9 3 8 2 9 4 0
5 8 0 0 9 6 5 7 7 8
1 9 1 4 0 5 7 4 2 5"""

# Convert string to matrix
matrix = [list(map(int, row.split())) for row in matrix_str.split('\n')]

# Rotate 3 times (equivalent to 630 degrees)
for _ in range(3):
    matrix = rotate_90_clockwise(matrix)

# Format output
result = '\n'.join(' '.join(map(str, row)) for row in matrix)
print(f"<<<{result}>>>")