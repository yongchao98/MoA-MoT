# Input matrix
matrix = [
    [5, 4],
    [0, 2],
    [9, 7],
    [6, 4],
    [9, 7],
    [1, 8],
    [2, 0],
    [5, 6],
    [1, 9],
    [5, 5]
]

# Perform max pooling
kernel_size = 2
result = []

# Process matrix in blocks of kernel_size
for i in range(0, len(matrix), kernel_size):
    block = matrix[i:i+kernel_size]  # Get 2x2 block
    # Flatten block and find maximum
    max_val = max(max(row) for row in block)
    result.append(max_val)

# Format output
output = "\n".join(f"{x:.2f}" for x in result)
print(output)