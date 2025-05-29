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
    
    result = []
    
    # Process each 2×2 block
    for i in range(0, rows - 1, kernel_size):
        # For each block, find maximum of all elements
        block_max = max(
            max(matrix[i][j], matrix[i+1][j])
            for j in range(cols)
        )
        result.append(f"{block_max:.2f}")
    
    # Print in required format
    print(f"<<<{chr(10).join(result)}>>>")

# Perform max pooling
max_pooling(matrix, 2)