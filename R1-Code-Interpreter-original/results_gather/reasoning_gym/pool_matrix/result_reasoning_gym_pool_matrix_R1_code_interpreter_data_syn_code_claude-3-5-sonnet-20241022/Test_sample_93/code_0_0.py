import numpy as np

# Input matrix
matrix = np.array([
    [2, 7, 3, 6, 0],
    [5, 6, 2, 5, 0],
    [4, 2, 9, 0, 7],
    [2, 7, 9, 5, 3],
    [5, 7, 9, 9, 1]
])

def max_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    output_rows = rows // kernel_size
    output_cols = cols // kernel_size
    output = np.zeros((output_rows, output_cols))
    
    for i in range(output_rows):
        for j in range(output_cols):
            # Extract the region and find maximum
            region = matrix[i*kernel_size:(i+1)*kernel_size, 
                          j*kernel_size:(j+1)*kernel_size]
            output[i,j] = np.max(region)
    
    return output

# Perform max pooling
result = max_pooling(matrix, 3)

# Format output
output_str = "<<<" + str(np.round(result, 2)) + ">>>"
print(output_str)