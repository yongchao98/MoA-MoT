# Input matrix
matrix = [
    [7, 5, 6, 8, 3],
    [7, 9, 3, 5, 9],
    [4, 1, 3, 1, 2]
]

# For kernel size 1, max pooling just returns the same matrix
# We just need to format the numbers
rows = len(matrix)
cols = len(matrix[0])

# Format output
result = ""
for i in range(rows):
    row = " ".join(f"{x:.2f}" for x in matrix[i])
    result += row
    if i < rows - 1:
        result += "\n"

print(f"<<<{result}>>>")