# Input matrix
matrix = [
    [8, 0, 3, 1, 8, 0, 6, 2],
    [1, 0, 7, 7, 6, 3, 9, 6],
    [4, 4, 0, 0, 0, 4, 0, 1],
    [3, 9, 6, 3, 9, 6, 0, 9],
    [1, 8, 0, 7, 8, 1, 6, 1],
    [9, 7, 0, 7, 5, 0, 9, 4],
    [9, 3, 6, 7, 1, 4, 6, 1],
    [3, 4, 4, 8, 4, 5, 2, 5],
    [9, 4, 4, 5, 0, 0, 0, 5]
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
        # Calculate average for each 2x2 block
        block_sum = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                block_sum += matrix[i*kernel_size + ki][j*kernel_size + kj]
        average = block_sum / (kernel_size * kernel_size)
        row.append(f"{average:.2f}")
    result.append(row)

# Format output
output = "\n".join([" ".join(row) for row in result])
print(output)