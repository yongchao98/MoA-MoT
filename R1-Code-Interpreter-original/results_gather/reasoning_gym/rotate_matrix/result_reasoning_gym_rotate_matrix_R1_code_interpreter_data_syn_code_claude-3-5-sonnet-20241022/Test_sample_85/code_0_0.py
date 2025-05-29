def rotate_180(matrix):
    # Convert the input string to a matrix
    matrix = [list(map(int, row.split())) for row in matrix.strip().split('\n')]
    
    # For 180-degree rotation, we can reverse the matrix and then reverse each row
    n = len(matrix)
    
    # First reverse the entire matrix
    matrix = matrix[::-1]
    
    # Then reverse each row
    for i in range(n):
        matrix[i] = matrix[i][::-1]
    
    # Format the output
    result = '\n'.join(' '.join(map(str, row)) for row in matrix)
    print(result)

# Input matrix
input_matrix = """2 2 3 6 2
8 6 1 8 9
5 5 7 0 2
5 6 1 2 6
6 6 6 0 2"""

rotate_180(input_matrix)