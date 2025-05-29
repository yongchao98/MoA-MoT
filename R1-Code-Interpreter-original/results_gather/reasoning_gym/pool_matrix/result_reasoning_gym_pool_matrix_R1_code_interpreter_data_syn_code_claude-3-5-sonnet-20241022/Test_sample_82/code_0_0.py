import numpy as np

# Input matrix
matrix = np.array([
    [3, 6, 2],
    [2, 8, 9],
    [3, 8, 4],
    [4, 0, 3],
    [5, 7, 0]
])

kernel_size = 2

def average_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    result_rows = (rows + kernel_size - 1) // kernel_size
    result_cols = (cols + kernel_size - 1) // kernel_size
    result = np.zeros((result_rows, result_cols))
    
    for i in range(0, rows, kernel_size):
        for j in range(0, cols, kernel_size):
            # Get the current window
            window = matrix[i:min(i+kernel_size, rows), j:min(j+kernel_size, cols)]
            # Calculate average for the window
            result[i//kernel_size, j//kernel_size] = np.mean(window)
    
    return result

# Perform average pooling
result = average_pooling(matrix, kernel_size)

# Format output to 2 decimal places
formatted_result = np.round(result, 2)
print(formatted_result)