import numpy as np

# Input matrix
matrix = [
    [3, 2, 6, 7, 1, 4, 9, 0, 0, 2],
    [7, 2, 2, 1, 1, 4, 1, 1, 8, 9],
    [4, 9, 9, 6, 1, 7, 1, 4, 5, 2],
    [0, 2, 5, 1, 3, 4, 8, 1, 3, 9],
    [2, 1, 0, 0, 2, 0, 0, 8, 3, 9],
    [0, 5, 0, 4, 5, 5, 1, 0, 7, 3],
    [5, 4, 2, 2, 3, 8, 3, 2, 3, 6],
    [4, 6, 0, 6, 3, 8, 7, 2, 6, 1]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Initialize output matrix
output = []

# Perform max pooling
for i in range(out_rows):
    row = []
    for j in range(out_cols):
        # Extract 3x3 region
        region = []
        for k in range(kernel_size):
            for l in range(kernel_size):
                region.append(matrix[i*kernel_size + k][j*kernel_size + l])
        # Find max value in region
        row.append(max(region))
    output.append(row)

# Format output with 2 decimal places
formatted_output = [[f"{x:.2f}" for x in row] for row in output]

# Print result in required format
print("<<<")
for row in formatted_output:
    print(" ".join(row))
print(">>>")