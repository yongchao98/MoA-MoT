# Input matrix
matrix = [
    [9, 5, 1, 1, 5, 2, 2, 6, 1],
    [2, 9, 0, 9, 8, 1, 9, 5, 6],
    [1, 1, 3, 3, 7, 2, 9, 4, 9],
    [4, 5, 3, 7, 1, 1, 0, 2, 9],
    [1, 1, 1, 4, 0, 3, 9, 3, 9],
    [6, 2, 5, 7, 1, 7, 6, 2, 8],
    [2, 8, 1, 4, 5, 9, 1, 8, 0],
    [4, 8, 0, 6, 8, 0, 0, 9, 7],
    [0, 0, 8, 6, 2, 8, 4, 8, 2]
]

# Since we're rotating 360 degrees (4 times 90 degrees),
# the matrix will return to its original state
# Let's format it properly for output

result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

# Remove the last newline and print
print(result.strip())