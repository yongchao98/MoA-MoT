# Input matrix
matrix = [
    [7, 9],
    [4, 4],
    [3, 3],
    [5, 2],
    [2, 5],
    [0, 3],
    [1, 5],
    [2, 6],
    [9, 4]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols

# Initialize output matrix
output = []

# Perform max pooling
for i in range(0, rows, kernel_size):
    row = []
    for j in range(cols):
        # Get the maximum value in the current kernel window
        max_val = max(matrix[i+k][j] for k in range(kernel_size))
        row.append(max_val)
    output.append(row)

# Format output
formatted_output = '\n'.join([' '.join(f'{x:.2f}' for x in row) for row in output])
print(formatted_output)