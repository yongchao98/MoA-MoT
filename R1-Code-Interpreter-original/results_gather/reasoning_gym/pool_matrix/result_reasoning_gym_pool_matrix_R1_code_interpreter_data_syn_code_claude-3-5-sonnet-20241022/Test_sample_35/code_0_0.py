import numpy as np

# Input matrix
matrix = np.array([
    [3, 5, 5],
    [7, 9, 9],
    [3, 6, 9]
])

def max_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = cols // kernel_size
    
    result = np.zeros((out_rows, out_cols))
    
    # Perform max pooling
    for i in range(out_rows):
        for j in range(out_cols):
            # Extract the region and find maximum
            region = matrix[i*kernel_size:(i+1)*kernel_size, 
                          j*kernel_size:(j+1)*kernel_size]
            result[i,j] = np.max(region)
    
    return result

# Perform max pooling with kernel size 2
result = max_pooling(matrix, 2)

# Format output to 2 decimal places
output = np.array2string(result, formatter={'float_kind':lambda x: "%.2f" % x})
print(output)