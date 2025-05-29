# Input matrix
matrix = [
    [7, 4, 4, 5, 8, 0, 9],
    [7, 5, 3, 0, 6, 4, 6],
    [4, 5, 7, 3, 1, 1, 9],
    [6, 4, 5, 8, 1, 8, 1],
    [6, 3, 4, 3, 9, 0, 8],
    [9, 5, 3, 2, 4, 4, 5]
]

# For kernel size 1, average pooling just returns the same matrix
# Format output to 2 decimal places
result = ""
for row in matrix:
    result += " ".join(f"{x:.2f}" for x in row) + "\n"

print(result.strip())