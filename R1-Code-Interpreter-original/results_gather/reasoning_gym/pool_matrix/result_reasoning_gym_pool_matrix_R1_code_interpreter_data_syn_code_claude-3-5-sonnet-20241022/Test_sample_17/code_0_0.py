# Input matrix
matrix = [
    [8, 5, 6, 4, 2, 5, 5],
    [1, 2, 4, 0, 2, 3, 6],
    [1, 6, 4, 8, 3, 8, 9]
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
output = " ".join(str(x) for x in result[0])
print(f"<<<{output}>>>")