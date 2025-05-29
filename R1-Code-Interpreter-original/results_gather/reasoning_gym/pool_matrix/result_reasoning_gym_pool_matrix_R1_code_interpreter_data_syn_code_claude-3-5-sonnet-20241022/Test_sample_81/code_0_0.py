# Input matrix
matrix = [
    [5, 2, 2, 6, 7, 0],
    [9, 0, 2, 7, 4, 5]
]

# For kernel size 1, max pooling just returns the same matrix
# We just need to format the numbers
rows = len(matrix)
cols = len(matrix[0])

# Format the output
result = ""
for i in range(rows):
    row = " ".join(f"{x:.2f}" for x in matrix[i])
    result += row
    if i < rows - 1:
        result += "\n"

print(f"<<<{result}>>>")