# Input matrix
input_str = """8 8 3 7 7 0 4 2 5
2 2 2 1 0 8 4 0 9
6 2 4 1 5 3 4 4 3
7 1 1 7 7 0 2 9 9"""

# Convert input to 2D list
matrix = [[int(x) for x in row.split()] for row in input_str.split('\n')]

def average_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = (cols + kernel_size - 1) // kernel_size
    
    result = []
    
    # Process each output cell
    for i in range(out_rows):
        row = []
        for j in range(out_cols):
            # Calculate pool boundaries
            r_start = i * kernel_size
            r_end = min(r_start + kernel_size, rows)
            c_start = j * kernel_size
            c_end = min(c_start + kernel_size, cols)
            
            # Calculate average of pool
            pool_sum = 0
            count = 0
            for r in range(r_start, r_end):
                for c in range(c_start, c_end):
                    pool_sum += matrix[r][c]
                    count += 1
            
            # Add average to result
            row.append(pool_sum / count)
        result.append(row)
    
    return result

# Get result and format output
result = average_pooling(matrix, 2)
output = '\n'.join([' '.join([f'{x:.2f}' for x in row]) for row in result])
print(f"<<<{output}>>>")