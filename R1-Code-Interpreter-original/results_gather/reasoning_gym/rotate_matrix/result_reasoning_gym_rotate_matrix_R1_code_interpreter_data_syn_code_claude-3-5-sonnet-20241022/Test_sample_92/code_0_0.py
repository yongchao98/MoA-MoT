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
input_str = """8 0 2 0 6 4
5 0 5 2 6 4
6 0 8 7 9 6
0 3 3 1 6 5
3 2 6 3 5 8
4 6 2 1 7 8"""

# Convert input to matrix
matrix = [list(map(int, row.split())) for row in input_str.splitlines()]

# Calculate effective rotation (810 % 360 = 90 degrees)
rotated_matrix = rotate_90_clockwise(matrix)

# Format output
output = '\n'.join(' '.join(map(str, row)) for row in rotated_matrix)
print(f"<<<{output}>>>")