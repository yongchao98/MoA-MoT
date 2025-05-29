# Input matrix
matrix = [
    [9, 1, 0, 4, 8, 9, 5],
    [4, 2, 5, 2, 2, 1, 6],
    [5, 2, 6, 9, 6, 7, 7],
    [3, 2, 9, 4, 8, 4, 0],
    [5, 9, 1, 5, 2, 6, 8],
    [8, 9, 7, 9, 8, 6, 9],
    [4, 8, 1, 9, 6, 5, 4]
]

# Since it's a 360-degree rotation, we just need to format the original matrix
result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

# Remove the last newline and print
print(result.strip())