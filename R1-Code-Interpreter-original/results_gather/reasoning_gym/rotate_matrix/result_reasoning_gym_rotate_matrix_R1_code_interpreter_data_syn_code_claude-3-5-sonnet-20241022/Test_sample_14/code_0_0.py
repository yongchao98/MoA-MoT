# Input matrix
matrix = [
    [6, 0, 9],
    [9, 3, 8],
    [2, 9, 8]
]

# For 0-degree rotation, we just format the output
result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

# Remove the last newline and print
print(result.strip())