# Input matrix
matrix = [
    [0, 0, 1, 4, 6, 5],
    [6, 6, 5, 2, 7, 4],
    [3, 6, 4, 1, 5, 9],
    [0, 6, 3, 2, 9, 3],
    [3, 3, 2, 0, 3, 2],
    [0, 3, 0, 0, 7, 3],
    [8, 7, 4, 4, 0, 0],
    [3, 3, 1, 9, 4, 5],
    [7, 0, 3, 5, 6, 1]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform average pooling
result = []
for i in range(0, rows - kernel_size + 1, kernel_size):
    row = []
    for j in range(0, cols - kernel_size + 1, kernel_size):
        # Calculate average of 2x2 region
        sum_region = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                sum_region += matrix[i + ki][j + kj]
        avg = sum_region / (kernel_size * kernel_size)
        row.append(round(avg, 2))
    if row:
        result.append(row)

# Print result in required format
output = "\n".join([" ".join(map(str, row)) for row in result])
print(f"<<<{output}>>>")