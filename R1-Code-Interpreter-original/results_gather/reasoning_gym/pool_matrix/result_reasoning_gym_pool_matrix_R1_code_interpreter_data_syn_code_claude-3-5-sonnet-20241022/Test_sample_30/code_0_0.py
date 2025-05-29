# Input matrix
matrix = [
    [2, 8, 1],
    [9, 3, 5],
    [1, 4, 5],
    [1, 9, 8]
]

# For kernel size 1, the output is the same as input
# Just format to 2 decimal places
rows = len(matrix)
cols = len(matrix[0])

# Format output
result = []
for i in range(rows):
    row = []
    for j in range(cols):
        row.append(f"{matrix[i][j]:.2f}")
    result.append(" ".join(row))

print("\n".join(result))