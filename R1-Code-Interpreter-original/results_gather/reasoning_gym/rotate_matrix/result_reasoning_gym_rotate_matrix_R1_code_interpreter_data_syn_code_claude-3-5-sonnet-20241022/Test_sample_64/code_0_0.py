def rotate_matrix(matrix):
    n = len(matrix)
    
    # Convert string input to 2D list
    matrix = [list(map(int, row.split())) for row in matrix.strip().split('\n')]
    
    # For 810 degrees = 90 degrees (810 % 360 = 90)
    # Only need one rotation
    
    # Step 1: Transpose
    for i in range(n):
        for j in range(i, n):
            matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]
    
    # Step 2: Reverse each row
    for i in range(n):
        matrix[i].reverse()
    
    # Format output
    result = []
    for row in matrix:
        result.append(' '.join(map(str, row)))
    return '\n'.join(result)

# Input matrix
input_matrix = """5 6 3 1 5 0 0
2 6 5 1 4 5 5
4 2 8 0 9 7 8
7 7 5 3 5 8 8
5 0 2 0 1 7 8
3 5 1 2 9 1 8
4 9 4 4 3 6 6"""

print(rotate_matrix(input_matrix))