# Input matrix
matrix = [
    [0, 3, 9, 0, 0],
    [0, 6, 2, 8, 9],
    [5, 0, 1, 1, 0],
    [8, 7, 9, 7, 1],
    [9, 2, 9, 1, 5],
    [0, 5, 5, 2, 2],
    [3, 6, 4, 0, 1],
    [7, 3, 1, 5, 5]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(0, rows, kernel_size):
    if i + kernel_size <= rows:  # Check if we have a complete kernel
        row = []
        for j in range(0, cols, kernel_size):
            if j + kernel_size <= cols:  # Check if we have a complete kernel
                # Get maximum value in current 2x2 region
                max_val = max(
                    matrix[i][j],
                    matrix[i][j+1],
                    matrix[i+1][j],
                    matrix[i+1][j+1]
                )
                row.append(max_val)
        if row:  # Only append if we have values in the row
            result.append(row)

# Format and print result
output = ""
for row in result:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"
print(output.strip())