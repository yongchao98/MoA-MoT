# Input matrix
matrix = [
    [5, 2],
    [8, 1],
    [9, 5],
    [7, 6],
    [4, 0],
    [6, 6],
    [6, 1]
]

# For kernel size 1, max pooling keeps the original values
# Just format the output with 2 decimal places
result = ""
for row in matrix:
    result += " ".join(f"{x:.2f}" for x in row) + "\n"

# Remove trailing newline and print
print(result.strip())