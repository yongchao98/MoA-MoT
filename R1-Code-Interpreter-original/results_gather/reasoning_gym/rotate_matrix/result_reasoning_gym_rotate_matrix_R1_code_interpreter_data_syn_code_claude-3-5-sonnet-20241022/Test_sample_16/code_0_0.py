# Input matrix
matrix = [
    [5, 2, 8, 8, 4, 8, 3, 6, 2, 2],
    [1, 7, 3, 3, 8, 8, 8, 2, 7, 3],
    [0, 1, 6, 1, 2, 7, 2, 2, 0, 8],
    [4, 5, 5, 1, 3, 4, 4, 0, 9, 4],
    [2, 5, 6, 2, 3, 6, 9, 4, 0, 5],
    [1, 4, 9, 4, 9, 7, 1, 1, 3, 2],
    [5, 7, 8, 5, 6, 1, 0, 5, 6, 6],
    [5, 0, 9, 6, 7, 6, 9, 4, 2, 5],
    [3, 7, 9, 3, 9, 2, 3, 2, 2, 4],
    [9, 6, 9, 8, 8, 6, 9, 6, 0, 7]
]

# Since 720 degrees is two full rotations, the matrix remains the same
# Let's format it properly for output
result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

print(result.strip())