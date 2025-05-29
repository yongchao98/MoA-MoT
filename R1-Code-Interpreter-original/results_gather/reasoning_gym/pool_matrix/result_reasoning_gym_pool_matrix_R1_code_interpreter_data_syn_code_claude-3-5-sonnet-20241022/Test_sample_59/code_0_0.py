import numpy as np

# Input matrix
matrix = np.array([
    [2, 2, 4, 4],
    [5, 2, 5, 1],
    [1, 8, 1, 2]
])

kernel_size = 3

def average_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    output_rows = rows // kernel_size
    output_cols = cols // kernel_size
    result = np.zeros((output_rows, output_cols))
    
    for i in range(output_rows):
        for j in range(output_cols):
            # Extract the region
            region = matrix[i*kernel_size:(i+1)*kernel_size, 
                          j*kernel_size:(j+1)*kernel_size]
            # Calculate average if the region is complete
            if region.shape == (kernel_size, kernel_size):
                result[i,j] = np.mean(region)
    
    return result

result = average_pooling(matrix, kernel_size)

# Format output to 2 decimal places
output = np.round(result, decimals=2)
print(output)