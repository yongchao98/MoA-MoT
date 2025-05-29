# Input matrix
matrix = [
    [7, 8, 8, 2, 1, 0, 6, 1, 3, 9],
    [0, 8, 0, 8, 1, 0, 5, 0, 2, 1],
    [6, 2, 6, 7, 2, 1, 6, 6, 3, 0],
    [0, 1, 0, 8, 0, 7, 8, 0, 5, 8],
    [2, 3, 4, 9, 6, 0, 3, 2, 3, 0],
    [2, 0, 0, 5, 7, 3, 3, 3, 3, 7],
    [6, 8, 9, 1, 6, 5, 4, 8, 3, 7]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Perform average pooling
result = []
for i in range(output_rows):
    row = []
    for j in range(output_cols):
        # Calculate average for each 3x3 region
        sum_region = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                sum_region += matrix[i*kernel_size + ki][j*kernel_size + kj]
        avg = sum_region / (kernel_size * kernel_size)
        row.append(round(avg, 2))
    result.append(row)

# Format output
output = "\n".join(" ".join(f"{x:.2f}" for x in row) for row in result)
print(f"<<<{output}>>>")