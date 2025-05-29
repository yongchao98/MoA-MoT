# Input matrix
matrix = [
    [3, 0, 5, 3, 4, 9, 8, 1],
    [5, 9, 6, 8, 7, 8, 5, 2],
    [8, 1, 6, 4, 8, 4, 6, 4],
    [9, 1, 6, 8, 8, 1, 6, 0],
    [4, 9, 8, 0, 9, 2, 6, 7],
    [0, 0, 2, 9, 2, 6, 0, 9]
]

# Kernel size
k = 2

# Get dimensions
rows = len(matrix)
cols = len(matrix[0])

# Initialize output matrix
output = []

# Process each 2x2 block
for i in range(0, rows, k):
    row = []
    for j in range(0, cols, k):
        # Get maximum value in current 2x2 block
        block_max = max(
            matrix[i][j],
            matrix[i][j+1],
            matrix[i+1][j],
            matrix[i+1][j+1]
        )
        row.append(block_max)
    output.append(row)

# Format output
formatted_output = '\n'.join([' '.join(f'{x:.2f}' for x in row) for row in output])
print(formatted_output)