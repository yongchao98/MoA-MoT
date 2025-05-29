import numpy as np

# Input matrix
matrix = [
    [2, 3, 5, 2, 7, 2, 4, 9],
    [5, 1, 0, 0, 7, 5, 9, 6],
    [5, 0, 6, 3, 1, 3, 6, 0],
    [7, 4, 7, 8, 8, 6, 7, 0],
    [4, 3, 7, 8, 4, 7, 5, 8]
]

def max_pooling(matrix, kernel_size):
    matrix = np.array(matrix)
    rows, cols = matrix.shape
    
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = cols // kernel_size
    
    result = np.zeros((out_rows, out_cols))
    
    # Perform max pooling
    for i in range(out_rows):
        for j in range(out_cols):
            # Extract the current region
            region = matrix[i*kernel_size:(i+1)*kernel_size, 
                          j*kernel_size:(j+1)*kernel_size]
            # Get maximum value in the region
            result[i, j] = np.max(region)
    
    return result

# Perform max pooling and format output
result = max_pooling(matrix, 3)
formatted_result = '\n'.join([' '.join([f'{x:.2f}' for x in row]) for row in result])
print(formatted_result)