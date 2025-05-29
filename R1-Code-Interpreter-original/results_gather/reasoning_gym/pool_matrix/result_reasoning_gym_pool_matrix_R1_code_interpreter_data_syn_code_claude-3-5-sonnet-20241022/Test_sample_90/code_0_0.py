# Input matrix
matrix = [
    [6, 0, 2, 6, 5],
    [9, 6, 3, 4, 3],
    [3, 0, 4, 1, 6],
    [6, 0, 5, 8, 7],
    [0, 8, 9, 2, 2],
    [5, 8, 3, 9, 5],
    [9, 9, 4, 2, 0],
    [7, 0, 4, 0, 5],
    [0, 2, 0, 2, 7],
    [1, 5, 7, 7, 6]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(out_rows):
    row = []
    for j in range(out_cols):
        # Extract the current kernel region
        values = []
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                r = i * kernel_size + ki
                c = j * kernel_size + kj
                if r < rows and c < cols:
                    values.append(matrix[r][c])
        # Find maximum value in the region
        max_val = max(values)
        row.append(f"{max_val:.2f}")
    result.append(row)

# Print result in required format
print("<<<")
for row in result:
    print(" ".join(row))
print(">>>")