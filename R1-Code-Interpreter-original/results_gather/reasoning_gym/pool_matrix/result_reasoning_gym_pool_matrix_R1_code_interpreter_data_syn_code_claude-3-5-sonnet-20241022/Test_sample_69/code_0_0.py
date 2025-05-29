# Input matrix
matrix = [
    [1, 6, 7, 9, 8, 4, 8, 5, 0, 5],
    [8, 1, 3, 8, 3, 3, 2, 8, 9, 3],
    [7, 0, 9, 7, 9, 8, 4, 3, 3, 7],
    [4, 8, 7, 6, 9, 4, 2, 7, 7, 7],
    [0, 4, 1, 8, 3, 1, 8, 4, 3, 2]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(0, rows, kernel_size):
    row = []
    for j in range(0, cols, kernel_size):
        # Get maximum value in current kernel window
        window_max = max(
            max(matrix[i+di][j+dj] 
                for dj in range(kernel_size))
            for di in range(kernel_size)
        )
        row.append(f"{window_max:.2f}")
    if row:  # Only append non-empty rows
        result.append(row)

# Format output
output = "\n".join([" ".join(row) for row in result])
print(f"<<<\n{output}\n>>>")