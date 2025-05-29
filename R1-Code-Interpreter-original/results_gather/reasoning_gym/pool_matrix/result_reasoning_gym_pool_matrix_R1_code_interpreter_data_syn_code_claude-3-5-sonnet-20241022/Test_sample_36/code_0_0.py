# Input matrix
matrix = [
    [2, 7, 4, 0, 8, 7, 4, 6, 3, 3],
    [8, 5, 1, 9, 8, 6, 1, 4, 5, 0],
    [7, 3, 5, 4, 7, 3, 7, 2, 5, 2],
    [9, 8, 0, 9, 0, 5, 7, 2, 9, 4],
    [4, 2, 2, 7, 7, 6, 0, 9, 1, 1],
    [0, 4, 3, 4, 7, 2, 5, 2, 6, 3],
    [4, 3, 4, 3, 8, 9, 2, 9, 6, 9],
    [4, 4, 9, 8, 7, 3, 2, 7, 2, 6],
    [3, 6, 5, 4, 8, 8, 7, 0, 5, 4]
]

# For kernel size 1, average pooling just returns the same matrix
# Format output to 2 decimal places
output = ""
for row in matrix:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"

print(output.strip())