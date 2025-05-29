# Input matrix
matrix = [
    [6, 3, 2, 3, 6, 7, 7, 3],
    [1, 5, 2, 4, 0, 9, 9, 4],
    [3, 2, 3, 9, 8, 7, 5, 0],
    [8, 3, 7, 3, 2, 0, 2, 2],
    [9, 8, 6, 8, 1, 8, 1, 6],
    [1, 6, 7, 9, 0, 9, 4, 8]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(output_rows):
    row = []
    for j in range(output_cols):
        # Get the values in the current 3x3 region
        values = []
        for x in range(kernel_size):
            for y in range(kernel_size):
                values.append(matrix[i*kernel_size + x][j*kernel_size + y])
        row.append(max(values))
    result.append(row)

# Format and print the output
output = ""
for row in result:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"
print(output.strip())