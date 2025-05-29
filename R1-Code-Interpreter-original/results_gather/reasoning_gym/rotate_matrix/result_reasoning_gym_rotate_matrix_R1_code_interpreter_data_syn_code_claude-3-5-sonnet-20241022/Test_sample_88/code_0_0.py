# Convert input string to matrix and perform rotation
def rotate_matrix_270(input_str):
    # Convert input string to matrix
    matrix = []
    rows = input_str.strip().split('\n')
    for row in rows:
        matrix.append(row.split())
    
    n = len(matrix)
    
    # For 270 degrees clockwise (or 90 degrees counterclockwise)
    # First transpose, then reverse each row
    # Create a new matrix for the result
    result = [[0 for _ in range(n)] for _ in range(n)]
    
    # Perform 270-degree clockwise rotation
    for i in range(n):
        for j in range(n):
            result[n-1-j][i] = matrix[i][j]
    
    # Format output
    output = []
    for row in result:
        output.append(' '.join(row))
    return '\n'.join(output)

# Input matrix
input_matrix = """3 3 9 7 1 5 3 3 4 5
5 2 1 5 4 2 4 4 4 3
3 8 8 1 6 3 3 0 1 0
7 6 6 9 6 6 8 6 1 0
0 6 5 4 9 3 2 6 6 9
5 7 4 5 6 4 4 6 4 8
8 5 0 5 8 3 8 3 5 5
8 4 2 8 9 1 0 5 3 9
1 3 2 6 0 9 2 7 1 7
0 1 9 3 2 8 5 5 3 7"""

print(rotate_matrix_270(input_matrix))