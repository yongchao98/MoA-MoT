# Input matrix
matrix = [
    [1, 7, 9, 9, 4],
    [5, 7, 9, 4, 4],
    [7, 4, 5, 1, 8],
    [2, 0, 6, 5, 7],
    [1, 2, 8, 8, 2]
]

# For 0-degree rotation, we just need to format the output
result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

# Remove the last newline and print
print(result.strip())