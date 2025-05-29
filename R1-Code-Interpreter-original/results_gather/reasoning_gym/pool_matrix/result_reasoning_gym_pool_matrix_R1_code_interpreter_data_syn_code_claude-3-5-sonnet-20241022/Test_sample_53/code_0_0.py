# Input matrix
matrix = [
    [0, 4, 3],
    [4, 3, 5],
    [1, 8, 3],
    [7, 6, 6],
    [0, 7, 9]
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
        # Calculate average for current window
        window_sum = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                window_sum += matrix[i + ki][j + kj]
        average = window_sum / (kernel_size * kernel_size)
        row.append(f"{average:.2f}")
    if row:
        result.append(row)

# Format and print result
output = "\n".join([" ".join(row) for row in result])
print(f"<<<{output}>>>")