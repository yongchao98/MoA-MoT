# Input matrix
matrix = [
    [1, 9, 6, 5, 4, 3, 5],
    [2, 9, 8, 6, 5, 6, 4],
    [5, 0, 5, 4, 3, 3, 5],
    [5, 8, 8, 4, 8, 4, 2],
    [5, 2, 2, 9, 6, 8, 2],
    [6, 0, 5, 4, 0, 4, 4],
    [9, 7, 4, 9, 7, 3, 6],
    [7, 9, 6, 8, 8, 9, 4],
    [4, 0, 1, 1, 5, 9, 3],
    [5, 0, 0, 2, 3, 4, 1]
]

kernel_size = 3
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
        # Calculate average for current 3x3 block
        block_sum = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                block_sum += matrix[i*kernel_size + ki][j*kernel_size + kj]
        avg = block_sum / (kernel_size * kernel_size)
        row.append(f"{avg:.2f}")
    result.append(row)

# Print result in required format
print("<<<")
for row in result:
    print(" ".join(row))
print(">>>")