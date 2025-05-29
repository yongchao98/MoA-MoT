# Input matrix
matrix = [
    [6, 6, 6, 3, 0],
    [9, 0, 8, 3, 3],
    [8, 3, 8, 7, 5],
    [6, 2, 9, 3, 8]
]

def average_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    result = []
    
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = cols // kernel_size
    
    for i in range(out_rows):
        row = []
        for j in range(out_cols):
            # Extract region and calculate average
            sum_val = 0
            for ki in range(kernel_size):
                for kj in range(kernel_size):
                    sum_val += matrix[i*kernel_size + ki][j*kernel_size + kj]
            avg = sum_val / (kernel_size * kernel_size)
            row.append(f"{avg:.2f}")
        result.append(row)
    
    # Print in required format
    print("<<<", end="")
    for i, row in enumerate(result):
        print(" ".join(row), end="")
    print(">>>")

# Perform average pooling with kernel size 3
average_pooling(matrix, 3)