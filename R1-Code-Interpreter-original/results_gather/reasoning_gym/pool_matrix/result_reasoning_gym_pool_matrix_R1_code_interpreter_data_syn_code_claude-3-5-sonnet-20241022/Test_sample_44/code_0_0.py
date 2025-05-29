# Input matrix
matrix = [
    [1, 7, 7, 0, 1, 2, 8, 7, 2],
    [8, 4, 5, 0, 3, 1, 4, 1, 2],
    [3, 6, 5, 2, 4, 5, 8, 6, 0],
    [4, 6, 4, 2, 3, 4, 8, 8, 8],
    [4, 9, 9, 9, 8, 5, 8, 6, 9],
    [4, 8, 7, 8, 7, 8, 3, 9, 8],
    [7, 2, 5, 5, 9, 0, 4, 4, 2]
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
        avg = sum_region / (kernel_size * kernel_size)
        row.append(round(avg, 2))
    result.append(row)

# Format output
output = "\n".join(" ".join(f"{x:.2f}" for x in row) for row in result)
print(f"<<<\n{output}\n>>>")