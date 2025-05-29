# Input matrix
matrix = [
    [8, 0, 1, 4, 5],
    [1, 0, 5, 0, 2],
    [1, 4, 7, 9, 7],
    [7, 6, 6, 3, 1],
    [0, 8, 6, 0, 0],
    [5, 7, 5, 3, 9],
    [5, 3, 7, 3, 9],
    [3, 8, 5, 0, 6]
]

def average_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Calculate output dimensions
    output_rows = rows // kernel_size
    output_cols = cols // kernel_size
    
    # Initialize output matrix
    output = []
    
    # Process each block
    for i in range(0, rows, kernel_size):
        if i + kernel_size <= rows:
            row = []
            for j in range(0, cols, kernel_size):
                if j + kernel_size <= cols:
                    # Calculate average of current block
                    block_sum = 0
                    for x in range(kernel_size):
                        for y in range(kernel_size):
                            block_sum += matrix[i + x][j + y]
                    average = block_sum / (kernel_size * kernel_size)
                    row.append(f"{average:.2f}")
            if row:
                output.append(row)
    
    # Print in required format
    print("<<<")
    for row in output:
        print(" ".join(row))
    print(">>>")

# Perform average pooling
average_pooling(matrix, 2)