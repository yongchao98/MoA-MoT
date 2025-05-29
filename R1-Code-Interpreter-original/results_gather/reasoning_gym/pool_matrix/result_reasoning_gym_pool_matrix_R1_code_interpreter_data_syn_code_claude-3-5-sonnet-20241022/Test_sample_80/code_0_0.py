# Input matrix
input_str = """4 8 1 4 5 7
8 4 9 0 4 8
6 0 3 5 6 5
9 8 3 1 5 2
4 0 8 0 0 8
7 4 2 8 2 1
3 4 2 3 7 0
7 8 2 7 9 1
7 8 9 6 5 9"""

# Convert input to matrix
matrix = [[int(num) for num in row.split()] for row in input_str.split('\n')]

# Parameters
kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])
output_rows = rows // kernel_size
output_cols = cols // kernel_size

# Perform max pooling
result = []
for i in range(0, rows, kernel_size):
    row = []
    for j in range(0, cols, kernel_size):
        # Get the current 3x3 block
        block = [matrix[x][y] for x in range(i, i + kernel_size) 
                for y in range(j, j + kernel_size)]
        # Find max in block
        max_val = max(block)
        row.append(f"{max_val:.2f}")
    if row:
        result.append(row)

# Format output
output = "\n".join([" ".join(row) for row in result])
print(f"<<<{output}>>>")