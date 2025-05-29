import numpy as np

# Input matrix
matrix = [
    [7, 6, 7],
    [8, 8, 3],
    [0, 7, 0],
    [7, 7, 1],
    [7, 2, 2],
    [1, 7, 4],
    [8, 9, 2]
]

def max_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Calculate output dimensions
    output_rows = rows // kernel_size
    output_cols = cols // kernel_size
    
    result = []
    
    # Process each block
    for i in range(0, output_rows * kernel_size, kernel_size):
        row = []
        for j in range(0, output_cols * kernel_size, kernel_size):
            # Get the block
            block = [matrix[x][y] for x in range(i, min(i + kernel_size, rows)) 
                    for y in range(j, min(j + kernel_size, cols))]
            # Find max in block
            max_val = max(block)
            row.append(f"{max_val:.2f}")
        if row:
            result.append(row)
    
    # Print in required format
    output = "\n".join([" ".join(row) for row in result])
    print(f"<<<{output}>>>")

# Perform max pooling
max_pooling(matrix, 2)