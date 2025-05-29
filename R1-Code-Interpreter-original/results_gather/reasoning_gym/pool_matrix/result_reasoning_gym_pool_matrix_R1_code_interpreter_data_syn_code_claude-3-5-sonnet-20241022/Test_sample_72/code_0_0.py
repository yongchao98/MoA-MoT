import numpy as np

# Input matrix
matrix = np.array([
    [5, 1, 9, 6, 8],
    [8, 4, 0, 3, 7]
])

kernel_size = 3

def average_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    # Calculate output dimensions
    out_rows = rows // kernel_size
    out_cols = cols // kernel_size
    
    result = []
    for i in range(0, rows, kernel_size):
        row_result = []
        for j in range(0, cols, kernel_size):
            # Get the region for pooling
            region = matrix[i:min(i+kernel_size, rows), j:min(j+kernel_size, cols)]
            # Calculate average of available elements
            avg = np.mean(region)
            if not np.isnan(avg):  # Only add if we have valid elements
                row_result.append(round(avg, 2))
        if row_result:  # Only add non-empty rows
            result.append(row_result)
    
    return np.array(result)

result = average_pooling(matrix, kernel_size)
print(result)