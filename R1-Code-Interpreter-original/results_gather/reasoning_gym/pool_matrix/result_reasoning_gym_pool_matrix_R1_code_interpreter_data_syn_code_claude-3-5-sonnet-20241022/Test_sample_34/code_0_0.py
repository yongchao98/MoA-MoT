# Input matrix
matrix = [
    [4, 7, 3, 5, 3, 0],
    [7, 6, 0, 2, 2, 2],
    [6, 0, 8, 8, 9, 9]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])
result = []

# Process each 2x2 region
for i in range(0, rows, kernel_size):
    row_result = []
    for j in range(0, cols, kernel_size):
        # Calculate average for current region
        sum_region = 0
        count = 0
        for ki in range(kernel_size):
            if i + ki < rows:  # Check if within bounds
                for kj in range(kernel_size):
                    if j + kj < cols:  # Check if within bounds
                        sum_region += matrix[i + ki][j + kj]
                        count += 1
        avg = sum_region / count
        row_result.append(f"{avg:.2f}")
    if row_result:  # Add row only if it contains elements
        result.append(row_result)

# Print result in required format
print(" ".join(result[0]))