import numpy as np

# Input matrix
matrix = np.array([
    [2, 2, 6, 1, 3, 9],
    [6, 1, 0, 1, 9, 0]
])

def max_pooling(matrix, kernel_size):
    rows, cols = matrix.shape
    result_rows = rows // kernel_size
    result_cols = cols // kernel_size
    
    result = np.zeros((result_rows, result_cols))
    
    for i in range(result_rows):
        for j in range(result_cols):
            window = matrix[
                i*kernel_size:(i+1)*kernel_size,
                j*kernel_size:(j+1)*kernel_size
            ]
            result[i,j] = np.max(window)
    
    return result

# Perform max pooling
result = max_pooling(matrix, 1)

# Format output
output = ""
for row in result:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"
print(output.strip())