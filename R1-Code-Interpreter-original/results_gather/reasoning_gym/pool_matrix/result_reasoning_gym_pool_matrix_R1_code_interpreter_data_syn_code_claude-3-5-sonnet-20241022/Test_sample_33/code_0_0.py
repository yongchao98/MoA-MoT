# Input matrix
matrix = [
    [5, 8, 5, 7, 8, 6],
    [6, 1, 8, 5, 1, 8],
    [9, 1, 7, 6, 8, 3],
    [8, 4, 0, 2, 3, 1],
    [6, 5, 0, 6, 3, 6],
    [1, 7, 7, 2, 4, 5],
    [8, 7, 2, 9, 0, 4],
    [2, 0, 4, 4, 8, 4]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Perform average pooling
result = []
for i in range(0, rows, kernel_size):
    row = []
    for j in range(0, cols, kernel_size):
        # Calculate average of 2x2 region
        sum_region = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                sum_region += matrix[i + ki][j + kj]
        average = sum_region / (kernel_size * kernel_size)
        row.append(f"{average:.2f}")
    if row:
        result.append(row)

# Print result in required format
print("<<<")
for row in result:
    print(" ".join(row))
print(">>>")