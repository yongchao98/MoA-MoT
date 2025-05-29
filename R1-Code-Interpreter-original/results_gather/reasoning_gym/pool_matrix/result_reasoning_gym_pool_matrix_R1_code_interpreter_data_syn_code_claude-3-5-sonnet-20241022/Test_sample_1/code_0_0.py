# Input matrix
matrix = [
    [4, 7, 8, 2, 2],
    [9, 9, 3, 6, 3],
    [3, 1, 7, 6, 0],
    [0, 6, 6, 0, 4]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(output_rows):
    row = []
    for j in range(output_cols):
        # Get the current 2x2 region
        region = []
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                region.append(matrix[i*kernel_size + ki][j*kernel_size + kj])
        # Find max in region
        row.append(max(region))
    result.append(row)

# Format output
output = "\n".join([" ".join([f"{x:.2f}" for x in row]) for row in result])
print(f"<<<{output}>>>")