# Input matrix
matrix = [
    [7, 3, 4, 8, 2, 8, 1, 7, 1],
    [7, 2, 9, 4, 2, 6, 4, 3, 6],
    [0, 6, 3, 3, 0, 6, 2, 9, 4],
    [2, 0, 4, 5, 0, 4, 6, 6, 4]
]

kernel_size = 2
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
        # Get the 2x2 region
        values = []
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                values.append(matrix[i*kernel_size + ki][j*kernel_size + kj])
        # Calculate average and round to 2 decimal places
        avg = round(sum(values) / len(values), 2)
        row.append(avg)
    result.append(row)

# Format output
output = "\n".join(" ".join(f"{x:.2f}" for x in row) for row in result)
print(output)