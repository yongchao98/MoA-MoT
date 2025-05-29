def max_pool(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = cols // kernel_size
    
    result = []
    
    # Process each pooling region
    for i in range(0, rows, kernel_size):
        if i + kernel_size <= rows:
            row_result = []
            for j in range(0, cols, kernel_size):
                if j + kernel_size <= cols:
                    # Get the current pooling region
                    pool = [matrix[x][y] for x in range(i, i + kernel_size) 
                           for y in range(j, j + kernel_size)]
                    # Find max in the region
                    max_val = max(pool)
                    row_result.append(max_val)
            if row_result:
                result.append(row_result)
    
    return result

# Input matrix
input_str = """6 5 5 0 1
4 9 3 0 2
5 4 5 6 0
4 6 5 6 8
0 7 4 3 6
1 5 4 6 5"""

# Convert input to matrix
matrix = [[int(num) for num in row.split()] for row in input_str.split('\n')]

# Perform max pooling
result = max_pool(matrix, 2)

# Format output
output = '\n'.join([' '.join(f'{x:.2f}' for x in row) for row in result])
print(f'<<<\n{output}\n>>>')