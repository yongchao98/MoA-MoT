# Input matrix
matrix = [
    [6, 1, 3, 0, 4, 9, 6, 1, 8],
    [2, 9, 6, 0, 6, 2, 0, 6, 3]
]

# For kernel size 1, we just need to format each number to 2 decimal places
rows = len(matrix)
cols = len(matrix[0])

# Create output string
output = ""
for i in range(rows):
    row_str = " ".join([f"{x:.2f}" for x in matrix[i]])
    output += row_str + "\n"

# Remove trailing newline and print
print(output.strip())