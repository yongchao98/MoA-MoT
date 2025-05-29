# Input matrix
matrix = [
    [5, 6, 7],
    [8, 9, 6],
    [8, 4, 4],
    [3, 3, 8],
    [4, 4, 6]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])

# Initialize output matrix
output_rows = rows // kernel_size
output_cols = cols // kernel_size
output = []

# Perform max pooling
for i in range(0, output_rows):
    row = []
    for j in range(0, output_cols):
        # Get values in current kernel region
        values = []
        for ki in range(kernel_size):
            for kj in range(kernel_size):
                values.append(matrix[i*kernel_size + ki][j*kernel_size + kj])
        # Add maximum value to output
        row.append(max(values))
    output.append(row)

# Format output
result = "\n".join([" ".join([f"{x:.2f}" for x in row]) for row in output])
print(f"<<<{result}>>>")