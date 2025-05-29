# Input matrix
matrix = [
    [8, 3, 5, 5],
    [7, 5, 8, 0],
    [8, 5, 1, 0]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])

# Initialize output matrix
output_rows = rows // kernel_size
output_cols = cols // kernel_size
output = []

# Perform max pooling
for i in range(0, output_rows):
    row = []
    for j in range(0, output_cols):
        # Get the current 2x2 region
        region = []
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                region.append(matrix[i*kernel_size + ki][j*kernel_size + kj])
        # Find max in region
        row.append(max(region))
    output.append(row)

# Format output
output_str = ' '.join(f'{x:.2f}' for x in output[0])
print(f'<<<{output_str}>>>')