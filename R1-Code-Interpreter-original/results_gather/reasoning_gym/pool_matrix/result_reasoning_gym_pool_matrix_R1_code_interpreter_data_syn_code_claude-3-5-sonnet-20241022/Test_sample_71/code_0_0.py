# Input matrix
matrix = [
    [7, 9, 3, 9, 6],
    [1, 5, 7, 7, 4],
    [6, 1, 4, 8, 1],
    [4, 5, 0, 1, 6],
    [0, 3, 5, 3, 8],
    [7, 6, 7, 1, 7],
    [8, 3, 5, 1, 8],
    [0, 1, 4, 4, 8],
    [9, 5, 5, 4, 5]
]

kernel_size = 3

# Calculate output dimensions
rows = len(matrix)
cols = len(matrix[0])
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform average pooling
result = []
for i in range(out_rows):
    row = []
    for j in range(out_cols):
        # Calculate average for current 3x3 region
        sum_region = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                sum_region += matrix[i*kernel_size + ki][j*kernel_size + kj]
        avg = sum_region / (kernel_size * kernel_size)
        row.append(f"{avg:.2f}")
    result.append(row)

# Format output
output = "\n".join([" ".join(row) for row in result])
print(f"<<<{output}>>>")