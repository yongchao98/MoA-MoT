# Input matrix
matrix = [
    [9, 5, 5, 6, 4, 0, 8, 9, 5],
    [2, 9, 4, 3, 6, 7, 9, 9, 7],
    [3, 6, 0, 6, 2, 0, 3, 5, 0],
    [5, 5, 5, 0, 7, 9, 6, 6, 4],
    [8, 8, 1, 8, 2, 0, 2, 6, 1]
]

kernel_size = 3

# Calculate output dimensions
rows = len(matrix)
cols = len(matrix[0])
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(out_rows):
    row = []
    for j in range(out_cols):
        # Extract the current 3x3 region
        region = []
        for r in range(kernel_size):
            for c in range(kernel_size):
                region.append(matrix[i*kernel_size + r][j*kernel_size + c])
        # Find max in the region
        row.append(max(region))
    result.append(row)

# Format output
output = " ".join(f"{x:.2f}" for x in result[0])
print(output)