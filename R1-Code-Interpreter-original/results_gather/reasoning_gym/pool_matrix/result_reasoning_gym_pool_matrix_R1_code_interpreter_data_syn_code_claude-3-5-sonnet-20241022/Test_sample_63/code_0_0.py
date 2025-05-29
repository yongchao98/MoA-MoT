# Input matrix
matrix = [
    [6, 9, 3, 6, 5, 2, 4, 7, 2, 1],
    [0, 4, 7, 2, 8, 1, 4, 3, 7, 3],
    [4, 3, 4, 2, 5, 4, 1, 8, 4, 6],
    [2, 6, 1, 5, 3, 6, 3, 6, 1, 0],
    [7, 3, 4, 1, 3, 3, 2, 8, 0, 3],
    [7, 3, 1, 7, 9, 4, 8, 4, 4, 5],
    [4, 3, 2, 6, 1, 4, 7, 4, 8, 7],
    [2, 0, 1, 0, 1, 8, 4, 3, 8, 8],
    [3, 0, 3, 4, 0, 3, 7, 4, 4, 8],
    [3, 1, 9, 8, 9, 2, 3, 4, 3, 9]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
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
        row.append(round(avg, 2))
    result.append(row)

# Format and print result
output = ""
for row in result:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"
print(output.strip())