# Input matrix
matrix = [
    [1, 7, 0, 2, 4],
    [0, 9, 0, 7, 6],
    [7, 8, 4, 2, 9],
    [5, 2, 8, 1, 1],
    [5, 0, 5, 4, 1]
]

# For 0-degree rotation, we just format the output
result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

# Remove the last newline and print
print(result.strip())