import numpy as np

# Input matrix
matrix = [
    [4, 3, 1, 3, 0, 4, 3, 8, 7, 7],
    [6, 9, 3, 7, 3, 3, 6, 5, 4, 5],
    [9, 1, 8, 7, 4, 5, 3, 0, 4, 9],
    [2, 8, 8, 6, 2, 0, 3, 4, 8, 3],
    [2, 2, 1, 2, 2, 9, 8, 1, 8, 9],
    [4, 2, 4, 6, 7, 5, 5, 6, 2, 5],
    [1, 8, 9, 1, 8, 0, 9, 3, 5, 9],
    [5, 0, 8, 0, 4, 2, 9, 7, 6, 6]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Initialize output matrix
output = []

# Perform average pooling
for i in range(out_rows):
    row = []
    for j in range(out_cols):
        # Get the current 3x3 region
        region_sum = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                region_sum += matrix[i*kernel_size + ki][j*kernel_size + kj]
        # Calculate average and round to 2 decimal places
        average = round(region_sum / (kernel_size * kernel_size), 2)
        row.append(average)
    output.append(row)

# Print the result in the required format
print("<<<")
for row in output:
    print(" ".join(str(x) for x in row))
print(">>>")