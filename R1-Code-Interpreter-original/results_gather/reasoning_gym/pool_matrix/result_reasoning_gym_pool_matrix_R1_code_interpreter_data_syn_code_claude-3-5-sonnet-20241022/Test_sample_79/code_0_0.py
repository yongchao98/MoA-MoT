import numpy as np

# Input matrix
matrix = [
    [4, 4, 4, 6, 3, 9, 2],
    [6, 2, 8, 7, 9, 7, 1],
    [0, 8, 2, 5, 6, 8, 5],
    [6, 5, 8, 6, 2, 6, 4],
    [3, 7, 9, 3, 1, 6, 8]
]

def average_pooling(matrix, kernel_size):
    matrix = np.array(matrix)
    rows, cols = matrix.shape
    pooled_rows = rows // kernel_size
    pooled_cols = cols // kernel_size
    
    result = np.zeros((pooled_rows, pooled_cols))
    
    for i in range(pooled_rows):
        for j in range(pooled_cols):
            # Extract the current window
            window = matrix[i*kernel_size:(i+1)*kernel_size, 
                          j*kernel_size:(j+1)*kernel_size]
            # Calculate average
            result[i, j] = np.mean(window)
    
    return result

# Perform pooling and round to 2 decimal places
result = average_pooling(matrix, 2)
result = np.round(result, 2)

# Format output
output = ' '.join([' '.join(map(str, row)) for row in result])
print(output)