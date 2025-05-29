def rotate_matrix(matrix):
    # Convert string input to 2D list
    matrix = [list(map(int, row.split())) for row in matrix.strip().split('\n')]
    n = len(matrix)
    
    # For 810 degrees = 90 degrees (810 % 360 = 90)
    # Create a new matrix for the transposition
    transposed = [[0 for _ in range(n)] for _ in range(n)]
    
    # Step 1: Transpose
    for i in range(n):
        for j in range(n):
            transposed[i][j] = matrix[j][i]
    
    # Step 2: Reverse each row
    for row in transposed:
        row.reverse()
    
    # Format output
    result = []
    for row in transposed:
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