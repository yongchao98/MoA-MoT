# Input matrix
matrix = [
    [7, 8, 8, 6, 6, 0],
    [9, 6, 5, 9, 1, 1],
    [0, 9, 8, 1, 1, 6],
    [4, 1, 4, 2, 0, 8],
    [0, 2, 8, 0, 5, 3],
    [0, 3, 9, 8, 6, 6],
    [8, 3, 6, 7, 2, 3],
    [0, 6, 7, 2, 3, 7],
    [5, 5, 7, 3, 2, 1],
    [1, 7, 4, 4, 6, 1]
]

# For kernel size 1, the output is the same as input
# Format the output with 2 decimal places
output = ""
for row in matrix:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"

print(output.strip())