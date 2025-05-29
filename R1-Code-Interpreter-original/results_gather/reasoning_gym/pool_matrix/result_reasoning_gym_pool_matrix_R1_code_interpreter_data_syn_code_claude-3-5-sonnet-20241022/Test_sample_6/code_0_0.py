import numpy as np

# Input matrix
matrix = [
    [1, 5, 3, 4, 2, 1, 7, 9],
    [0, 1, 4, 3, 1, 7, 2, 1],
    [8, 2, 3, 3, 4, 4, 3, 4],
    [7, 5, 3, 7, 4, 9, 2, 3],
    [5, 9, 7, 1, 4, 5, 7, 4],
    [6, 2, 6, 5, 0, 5, 0, 5],
    [3, 3, 6, 5, 2, 5, 6, 7],
    [1, 1, 7, 4, 9, 3, 0, 3]
]

kernel_size = 3
matrix = np.array(matrix)
rows, cols = matrix.shape
output_rows = rows // kernel_size
output_cols = cols // kernel_size
output = np.zeros((output_rows, output_cols))

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current 3x3 region
        region = matrix[i*kernel_size:(i+1)*kernel_size, 
                       j*kernel_size:(j+1)*kernel_size]
        # Calculate average
        output[i,j] = np.mean(region)

# Format output to string with 2 decimal places
output_str = "<<< "
for row in output:
    output_str += " ".join([f"{x:.2f}" for x in row]) + "\n"
output_str += ">>>"
print(output_str)