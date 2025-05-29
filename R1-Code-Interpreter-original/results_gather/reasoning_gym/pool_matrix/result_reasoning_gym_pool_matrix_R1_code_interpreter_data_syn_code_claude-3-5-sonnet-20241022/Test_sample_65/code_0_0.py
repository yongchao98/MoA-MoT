# Input matrix
matrix = [
    [3, 9, 3, 5, 7, 3, 5, 4],
    [8, 2, 3, 4, 4, 4, 4, 3],
    [0, 8, 3, 3, 8, 7, 9, 9],
    [8, 3, 5, 9, 0, 3, 3, 4],
    [7, 1, 9, 5, 1, 1, 7, 7]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Perform average pooling
result = []
for i in range(output_rows):
    row = []
    for j in range(output_cols):
        # Calculate average for each 2x2 region
        sum_region = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                sum_region += matrix[i*kernel_size + ki][j*kernel_size + kj]
        average = sum_region / (kernel_size * kernel_size)
        row.append(round(average, 2))
    result.append(row)

# Format output
output = "\n".join(" ".join(f"{x:.2f}" for x in row) for row in result)
print(output)