# Input matrix
matrix = [
    [6, 9, 7],
    [4, 1, 6]
]

# For kernel size 1, average pooling just maintains the same values
# but we need to format to 2 decimal places
rows = len(matrix)
cols = len(matrix[0])

# Format output
result = ""
for i in range(rows):
    for j in range(cols):
        result += f"{matrix[i][j]:.2f}"
        if j < cols - 1:
            result += " "
    if i < rows - 1:
        result += "\n"

print(result)