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
for i in range(output_rows):
    # Process each 2x2 region
    region_values = []
    for ki in range(kernel_size):
        for kj in range(kernel_size):
            region_values.append(matrix[i*kernel_size + ki][kj])
    output.append(max(region_values))

# Format output
result = "\n".join([f"{x:.2f}" for x in output])
print(f"<<<{result}>>>")