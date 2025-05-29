# Input matrix
matrix = [
    [9, 9, 4],
    [9, 5, 9],
    [3, 3, 8],
    [5, 5, 2],
    [3, 3, 7],
    [5, 1, 4],
    [4, 2, 9],
    [6, 7, 6],
    [7, 9, 9],
    [3, 6, 9]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform average pooling
result = []
for i in range(0, rows, kernel_size):
    if i + kernel_size <= rows:  # Only process complete 2x2 blocks
        row = []
        for j in range(0, cols, kernel_size):
            if j + kernel_size <= cols:
                # Calculate average of current kernel window
                sum_val = 0
                for ki in range(kernel_size):
                    for kj in range(kernel_size):
                        sum_val += matrix[i + ki][j + kj]
                avg = sum_val / (kernel_size * kernel_size)
                row.append(f"{avg:.2f}")
        if row:  # Only append non-empty rows
            result.append(row)

# Print result in required format
print("<<<")
for row in result:
    print(" ".join(row))
print(">>>")