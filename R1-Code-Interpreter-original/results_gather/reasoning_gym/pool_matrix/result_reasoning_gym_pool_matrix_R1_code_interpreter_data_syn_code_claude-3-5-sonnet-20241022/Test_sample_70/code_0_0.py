# Input matrix
matrix = [
    [0, 3, 2, 2, 2, 9, 8],
    [8, 9, 7, 4, 7, 0, 6],
    [2, 7, 3, 0, 0, 0, 6],
    [7, 9, 7, 4, 9, 3, 2]
]

# For kernel size 1, average pooling just returns the same matrix
# We'll format the output to 2 decimal places
rows = len(matrix)
cols = len(matrix[0])

result = []
for i in range(rows):
    row = []
    for j in range(cols):
        # Format to 2 decimal places
        value = "{:.2f}".format(matrix[i][j])
        row.append(value)
    result.append(" ".join(row))

print("\n".join(result))