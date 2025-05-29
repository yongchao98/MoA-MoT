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
matrix_str = """9 8 0 3 4 3 2 5 4
7 0 8 2 0 3 9 7 0
2 9 5 1 4 5 2 8 8
1 6 0 5 5 3 2 4 2
8 4 4 7 2 9 3 7 6
8 1 2 9 3 9 9 0 1
7 9 4 1 3 3 7 7 1
7 5 9 7 4 8 9 3 1
0 2 9 6 2 7 1 6 3"""

# Convert string to matrix
matrix = [list(map(int, row.split())) for row in matrix_str.split('\n')]

# Calculate effective rotation (810 degrees = 90 degrees)
rotated_matrix = rotate_90_clockwise(matrix)

# Format output
output = '\n'.join(' '.join(map(str, row)) for row in rotated_matrix)
print(f"<<<{output}>>>")