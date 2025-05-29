# Input matrix
matrix = [
    [0, 9, 2, 9, 8, 0],
    [8, 5, 9, 2, 2, 3],
    [1, 1, 8, 3, 4, 1],
    [8, 5, 9, 1, 3, 6]
]

# For kernel size 1, we just need to format each number
rows = len(matrix)
cols = len(matrix[0])

# Create output string
output = ""
for i in range(rows):
    row_str = " ".join(f"{x:.2f}" for x in matrix[i])
    output += row_str + "\n"

# Remove trailing newline and print
print(output.strip())