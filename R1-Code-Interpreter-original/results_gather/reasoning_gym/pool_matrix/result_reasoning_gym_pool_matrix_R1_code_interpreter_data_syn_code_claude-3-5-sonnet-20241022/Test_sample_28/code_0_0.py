# Input matrix
matrix = [
    [1, 6, 3, 8, 9, 1],
    [2, 1, 7, 8, 5, 2],
    [8, 6, 3, 2, 5, 0],
    [3, 4, 0, 6, 6, 9],
    [5, 7, 0, 2, 5, 8],
    [7, 9, 4, 0, 7, 0]
]

kernel_size = 3
n = len(matrix)
output_size = n // kernel_size

# Initialize output matrix
output = [[0 for _ in range(output_size)] for _ in range(output_size)]

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Calculate average for each 3x3 region
        sum_region = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                sum_region += matrix[i*kernel_size + ki][j*kernel_size + kj]
        output[i][j] = round(sum_region / (kernel_size * kernel_size), 2)

# Format and print output
result = ""
for row in output:
    result += " ".join(f"{x:.2f}" for x in row) + "\n"
print(result.strip())