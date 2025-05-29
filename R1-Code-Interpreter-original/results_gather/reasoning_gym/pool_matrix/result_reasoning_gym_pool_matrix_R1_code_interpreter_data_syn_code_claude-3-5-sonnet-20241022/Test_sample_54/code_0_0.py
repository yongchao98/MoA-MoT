# Input matrix
matrix = [
    [4, 0, 3, 8, 9, 9, 9, 1, 7, 8],
    [9, 4, 4, 1, 0, 7, 0, 5, 8, 0]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])
result = []

# Process each 2x2 region
for j in range(0, cols, kernel_size):
    # Calculate average for current region
    sum_region = (matrix[0][j] + matrix[0][j+1] + 
                 matrix[1][j] + matrix[1][j+1])
    avg = sum_region / (kernel_size * kernel_size)
    result.append(f"{avg:.2f}")

# Format output
print(" ".join(result))