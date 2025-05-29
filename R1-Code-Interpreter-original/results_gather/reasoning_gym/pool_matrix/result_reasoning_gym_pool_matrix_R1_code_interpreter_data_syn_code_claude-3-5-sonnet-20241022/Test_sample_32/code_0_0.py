# Input matrix
matrix = [
    [4, 6, 7],
    [6, 6, 1],
    [4, 3, 5],
    [0, 8, 6],
    [8, 2, 0],
    [1, 4, 2],
    [6, 4, 5],
    [0, 5, 9],
    [6, 2, 6]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(out_rows):
    row_start = i * kernel_size
    max_val = float('-inf')
    # Find maximum in the 3x3 region
    for r in range(row_start, row_start + kernel_size):
        for c in range(cols):
            max_val = max(max_val, matrix[r][c])
    result.append([round(max_val, 2)])

# Print result in required format
output = "\n".join(" ".join(str(x) for x in row) for row in result)
print(f"<<<{output}>>>")