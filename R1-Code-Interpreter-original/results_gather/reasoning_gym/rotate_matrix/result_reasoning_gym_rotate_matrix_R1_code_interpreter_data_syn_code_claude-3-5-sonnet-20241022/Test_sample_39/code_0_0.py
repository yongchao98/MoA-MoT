# Input matrix
matrix = [
    [8, 9, 0, 4, 2, 3, 4],
    [1, 3, 8, 0, 9, 3, 8],
    [0, 8, 5, 7, 2, 5, 2],
    [0, 5, 9, 0, 6, 9, 3],
    [2, 5, 6, 7, 1, 0, 3],
    [9, 2, 2, 3, 8, 5, 4],
    [9, 0, 4, 5, 3, 6, 0]
]

# Step 1: Transpose the matrix
n = len(matrix)
for i in range(n):
    for j in range(i, n):
        matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]

# Step 2: Reverse each row
for i in range(n):
    matrix[i].reverse()

# Format output
output = ""
for row in matrix:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())