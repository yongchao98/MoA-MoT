def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix to store the rotated values
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Fill the result matrix
    # For clockwise rotation:
    # first row becomes last column
    # second row becomes second-to-last column
    # and so on...
    for i in range(n):
        for j in range(n):
            result[j][n-1-i] = matrix[i][j]
    
    return result

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

# Rotate matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format output
output = '\n'.join(' '.join(map(str, row)) for row in rotated_matrix)
print(f"<<<{output}>>>")