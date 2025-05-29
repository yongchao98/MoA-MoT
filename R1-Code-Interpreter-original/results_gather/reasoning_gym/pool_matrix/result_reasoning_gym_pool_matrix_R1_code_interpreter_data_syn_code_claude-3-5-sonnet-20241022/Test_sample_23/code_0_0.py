# Input matrix
matrix = [
    [7, 8, 5, 7, 5],
    [1, 6, 0, 3, 8],
    [7, 1, 5, 7, 8]
]

# For kernel size 1, max pooling just returns the same matrix
# We just need to format the output
rows = len(matrix)
cols = len(matrix[0])

# Format output
result = ""
for i in range(rows):
    row = " ".join(f"{x:.2f}" for x in matrix[i])
    result += row
    if i < rows - 1:
        result += "\n"

print(result)