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
        # Get the maximum value from the 2×2 block
        max_val = max(
            max(matrix[i][j:j+kernel_size]) for j in range(0, cols, kernel_size)
        )
        result.append(f"{max_val:.2f}")
    
    # Print in required format
    print(f"<<<{chr(10).join(result)}>>>")

# Perform max pooling
max_pooling(matrix, 2)