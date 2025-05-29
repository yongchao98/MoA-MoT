# Input matrix
matrix = [
    [3, 1, 8, 3, 1, 1, 1, 3, 9],
    [4, 3, 8, 4, 1, 0, 5, 4, 1],
    [2, 4, 1, 3, 3, 1, 0, 9, 1],
    [6, 3, 6, 9, 1, 5, 2, 6, 9],
    [1, 0, 0, 4, 5, 5, 6, 1, 3],
    [2, 9, 3, 8, 4, 4, 6, 9, 9],
    [5, 2, 6, 0, 1, 4, 4, 4, 2],
    [9, 2, 7, 8, 1, 0, 1, 2, 6],
    [5, 1, 9, 7, 0, 0, 5, 8, 6],
    [4, 8, 1, 1, 0, 1, 7, 6, 3]
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
    if i + kernel_size <= rows:  # Only process complete blocks
        row = []
        for j in range(0, cols, kernel_size):
            if j + kernel_size <= cols:  # Only process complete blocks
                # Get maximum value in current kernel window
                max_val = max(
                    matrix[i][j],
                    matrix[i][j+1],
                    matrix[i+1][j],
                    matrix[i+1][j+1]
                )
                row.append(max_val)
        if row:  # Only append non-empty rows
            result.append(row)

# Format and print the result
output = "<<< "
for row in result:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"
output += ">>>"
print(output)