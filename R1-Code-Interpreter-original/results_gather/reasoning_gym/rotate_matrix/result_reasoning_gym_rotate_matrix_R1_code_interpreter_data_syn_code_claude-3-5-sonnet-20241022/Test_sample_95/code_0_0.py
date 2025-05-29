# Input matrix
matrix = [
    [9, 9],
    [9, 1]
]

# For 0 degree rotation, we just return the matrix in required format
result = ""
for row in matrix:
    result += " ".join(map(str, row)) + "\n"

# Remove the last newline and print
print(result.strip())