import numpy as np

def max_pooling(matrix, kernel_size):
    rows, cols = len(matrix), len(matrix[0])
    output_rows = rows // kernel_size
    output_cols = cols // kernel_size
    output_matrix = np.zeros((output_rows, output_cols))
    
    for i in range(output_rows):
        for j in range(output_cols):
            sub_matrix = matrix[i*kernel_size:(i+1)*kernel_size, j*kernel_size:(j+1)*kernel_size]
            output_matrix[i, j] = np.max(sub_matrix)
    
    return output_matrix

input_matrix = np.array([
    [3, 0, 5, 3, 4, 9, 8, 1],
    [5, 9, 6, 8, 7, 8, 5, 2],
    [8, 1, 6, 4, 8, 4, 6, 4],
    [9, 1, 6, 8, 8, 1, 6, 0],
    [4, 9, 8, 0, 9, 2, 6, 7],
    [0, 0, 2, 9, 2, 6, 0, 9]
])

output_matrix = max_pooling(input_matrix, 2)
print(np.round(output_matrix, 2))