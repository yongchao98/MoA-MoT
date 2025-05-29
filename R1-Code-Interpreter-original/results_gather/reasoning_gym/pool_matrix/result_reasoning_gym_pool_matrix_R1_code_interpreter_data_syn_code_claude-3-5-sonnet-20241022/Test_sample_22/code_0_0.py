# Input matrix
matrix = [
    [5, 2, 7, 1, 0, 8, 7, 1, 9],
    [6, 9, 7, 4, 4, 4, 6, 4, 0],
    [6, 8, 9, 7, 3, 2, 1, 2, 6],
    [9, 8, 1, 7, 7, 2, 7, 7, 8]
]

def max_pooling(matrix, kernel_size):
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Calculate output dimensions
    output_rows = rows // kernel_size
    output_cols = cols // kernel_size
    
    result = []
    
    # Process each region
    for i in range(output_rows):
        row = []
        for j in range(output_cols):
            # Get the current region
            region = []
            for r in range(kernel_size):
                for c in range(kernel_size):
                    row_idx = i * kernel_size + r
                    col_idx = j * kernel_size + c
                    if row_idx < rows and col_idx < cols:
                        region.append(matrix[row_idx][col_idx])
            
            # Find maximum in the region
            max_val = max(region)
            row.append(f"{max_val:.2f}")
        result.append(row)
    
    # Format output
    output = "\n".join([" ".join(row) for row in result])
    print(f"<<<{output}>>>")

# Perform max pooling with kernel size 3
max_pooling(matrix, 3)