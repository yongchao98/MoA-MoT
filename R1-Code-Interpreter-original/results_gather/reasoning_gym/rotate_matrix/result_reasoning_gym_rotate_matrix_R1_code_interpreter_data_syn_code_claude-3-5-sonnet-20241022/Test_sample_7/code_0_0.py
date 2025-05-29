# Create the input matrix
matrix = [
    [4, 6, 0, 3, 3],
    [6, 9, 0, 4, 5],
    [2, 0, 4, 1, 5],
    [0, 1, 6, 5, 6],
    [1, 6, 8, 9, 4]
]

# For 360 degrees (4 rotations of 90 degrees), we get back the original matrix
# Let's format it properly for output
result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

# Remove the last newline and print
print(result.strip())