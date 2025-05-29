# Input matrix
input_str = """9 9 9 0 4 6
3 0 4 8 1 4
9 0 4 4 4 9
6 5 7 8 8 8
7 1 9 7 6 9"""

# Convert input to matrix
matrix = [[int(num) for num in row.split()] for row in input_str.split('\n')]

# Get dimensions
rows = len(matrix)
cols = len(matrix[0])
kernel_size = 2

# Initialize output matrix
output_rows = rows // kernel_size
output_cols = cols // kernel_size
output = []

# Perform average pooling
for i in range(0, rows, kernel_size):
    output_row = []
    for j in range(0, cols, kernel_size):
        # Calculate average of current 2x2 block
        block_sum = 0
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                if i + ki < rows and j + kj < cols:
                    block_sum += matrix[i + ki][j + kj]
        average = block_sum / (kernel_size * kernel_size)
        output_row.append(f"{average:.2f}")
    if output_row:
        output.append(output_row)

# Format output
result = "\n".join([" ".join(row) for row in output])
print(f"<<<{result}>>>")