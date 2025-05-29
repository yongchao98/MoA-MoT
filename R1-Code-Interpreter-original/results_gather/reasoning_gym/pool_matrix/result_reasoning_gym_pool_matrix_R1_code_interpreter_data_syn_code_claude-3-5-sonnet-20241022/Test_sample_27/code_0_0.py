import numpy as np

# Input matrix
matrix = np.array([
    [1, 7, 5],
    [9, 0, 3],
    [0, 3, 7],
    [9, 0, 9],
    [4, 5, 2]
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

# Get result and format
result = max_pooling(matrix, 3)
# Format output with 2 decimal places
formatted_result = np.array2string(result, formatter={'float_kind':lambda x: "%.2f" % x})
print(formatted_result)