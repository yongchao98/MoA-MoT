def max_pooling(matrix, kernel_size):
    # Convert input string to matrix
    matrix = [[int(x) for x in row.split()] for row in matrix.strip().split('\n')]
    
    rows = len(matrix)
    cols = len(matrix[0])
    
    # For kernel_size 1, output is same as input
    if kernel_size == 1:
        result = matrix
    
    # Format output with 2 decimal places
    formatted_result = []
    for row in result:
        formatted_row = [f"{x:.2f}" for x in row]
        formatted_result.append(' '.join(formatted_row))
    
    print('\n'.join(formatted_result))

# Input matrix
input_matrix = """7 0 8 4 1
7 1 3 8 9
4 9 9 6 7
9 1 9 4 5
7 2 2 2 0
0 7 4 2 2
4 8 5 6 4
1 9 2 0 5
7 5 7 8 4
2 1 6 9 5"""

max_pooling(input_matrix, 1)