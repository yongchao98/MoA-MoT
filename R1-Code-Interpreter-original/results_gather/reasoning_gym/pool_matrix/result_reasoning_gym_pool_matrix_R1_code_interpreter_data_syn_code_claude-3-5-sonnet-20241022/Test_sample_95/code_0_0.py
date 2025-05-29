# Input matrix
matrix = [
    [9, 2, 4],
    [3, 9, 4]
]

# For kernel size 1, average pooling just returns the same matrix
# but with numbers formatted to 2 decimal places
rows = len(matrix)
cols = len(matrix[0])

# Create output string
output = ""
for i in range(rows):
    row_str = " ".join(f"{x:.2f}" for x in matrix[i])
    output += row_str
    if i < rows - 1:
        output += "\n"

print(output)