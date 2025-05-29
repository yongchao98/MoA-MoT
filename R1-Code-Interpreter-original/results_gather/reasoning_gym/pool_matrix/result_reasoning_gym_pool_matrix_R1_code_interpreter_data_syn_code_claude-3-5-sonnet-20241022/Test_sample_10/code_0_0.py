import numpy as np

# Input matrix
matrix = np.array([
    [6, 3],
    [7, 4],
    [6, 9]
])

def max_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = cols // kernel_size
    
    result = np.zeros((out_rows, out_cols))
    
    # For this specific case, we only have one region
    for i in range(out_rows):
        for j in range(out_cols):
            # Extract the region and find maximum
            region = matrix[i*kernel_size:(i+1)*kernel_size, 
                          j*kernel_size:(j+1)*kernel_size]
            result[i,j] = np.max(region)
    
    return result

result = max_pooling(matrix, 3)
# Format to 2 decimal places
formatted_result = np.round(result, 2)
print(formatted_result)