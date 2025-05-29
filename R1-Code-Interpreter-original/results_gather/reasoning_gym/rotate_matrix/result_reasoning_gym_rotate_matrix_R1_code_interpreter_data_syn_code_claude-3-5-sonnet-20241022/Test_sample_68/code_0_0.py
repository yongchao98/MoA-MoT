def rotate_matrix_270(matrix):
    # Convert the input string to a matrix
    rows = matrix.strip().split('\n')
    matrix = [list(map(int, row.split())) for row in rows]
    
    n = len(matrix)
    # For 270 degrees clockwise, elements move according to the pattern:
    # (i,j) -> (j, n-1-i)
    # Create a new matrix
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            # For 270 degrees clockwise, new_row = j, new_col = n-1-i
            rotated[j][n-1-i] = matrix[i][j]
    
    # Format the output
    result = []
    for row in rotated:
        result.append(' '.join(map(str, row)))
    
    print('<<<' + '\n'.join(result) + '>>>')

# Input matrix
input_matrix = """7 7 2 6
5 6 8 1
8 1 1 7
4 0 3 5"""

rotate_matrix_270(input_matrix)