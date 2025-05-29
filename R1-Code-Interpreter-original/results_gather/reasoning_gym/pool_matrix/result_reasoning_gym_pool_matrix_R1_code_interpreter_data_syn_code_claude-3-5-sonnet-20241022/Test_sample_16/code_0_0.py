# Input matrix
matrix = [
    [0, 3, 1, 4, 0, 0, 6],
    [6, 2, 0, 6, 6, 0, 7],
    [8, 2, 2, 5, 2, 5, 9],
    [2, 6, 2, 1, 5, 6, 0],
    [5, 2, 6, 6, 5, 6, 2],
    [2, 5, 4, 7, 3, 2, 0],
    [3, 3, 6, 0, 2, 1, 6],
    [9, 3, 9, 1, 7, 4, 2],
    [1, 6, 3, 5, 0, 7, 4],
    [2, 1, 9, 3, 0, 3, 0]
]

kernel_size = 1

# For kernel_size = 1, output is same as input
# Format output with 2 decimal places
output = ""
for row in matrix:
    output += " ".join(f"{x:.2f}" for x in row) + "\n"

print(output.strip())