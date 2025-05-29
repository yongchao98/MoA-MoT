# Input matrix
matrix = [
    [0, 9, 6, 2, 9, 6, 0, 2],
    [2, 5, 9, 4, 0, 9, 8, 4],
    [3, 8, 0, 5, 6, 9, 7, 3],
    [6, 6, 9, 0, 1, 9, 4, 4],
    [5, 8, 4, 2, 9, 9, 8, 1]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(out_rows):
    row = []
    for j in range(out_cols):
        # Get the 2x2 region
        region = [
            matrix[i*kernel_size + x][j*kernel_size + y]
            for x in range(kernel_size)
            for y in range(kernel_size)
        ]
        # Find max in region
        max_val = max(region)
        row.append(max_val)
    result.append(row)

# Format output
output = "\n".join(" ".join(f"{x:.2f}" for x in row) for row in result)
print(f"<<<{output}>>>")