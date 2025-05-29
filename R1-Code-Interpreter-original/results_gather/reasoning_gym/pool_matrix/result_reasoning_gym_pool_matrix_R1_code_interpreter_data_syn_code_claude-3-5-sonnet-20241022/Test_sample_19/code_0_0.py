import numpy as np

# Input matrix
matrix = [
    [9, 5, 0, 9, 5, 0, 5, 4, 9, 5],
    [5, 4, 4, 3, 6, 4, 6, 9, 1, 4],
    [7, 0, 8, 8, 6, 4, 5, 9, 2, 3],
    [7, 7, 3, 8, 2, 5, 7, 1, 3, 5],
    [8, 9, 6, 5, 7, 5, 1, 4, 1, 1]
]

def average_pooling(matrix, kernel_size):
    matrix = np.array(matrix)
    rows, cols = matrix.shape
    
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = cols // kernel_size
    
    result = np.zeros((out_rows, out_cols))
    
    # Perform pooling
    for i in range(out_rows):
        for j in range(out_cols):
            # Extract the current window
            window = matrix[i*kernel_size:(i+1)*kernel_size, 
                          j*kernel_size:(j+1)*kernel_size]
            # Calculate average
            result[i, j] = np.mean(window)
    
    return result

# Perform average pooling with kernel size 3
result = average_pooling(matrix, 3)

# Format output to 2 decimal places
formatted_result = np.round(result, 2)
print(formatted_result)