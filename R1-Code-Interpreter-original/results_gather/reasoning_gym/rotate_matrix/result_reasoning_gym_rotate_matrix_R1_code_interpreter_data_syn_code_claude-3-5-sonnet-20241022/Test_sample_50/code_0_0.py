# Input matrix
matrix = [
    [7, 8, 4, 2],
    [2, 4, 2, 0],
    [6, 7, 9, 1],
    [3, 4, 8, 2]
]

# Get dimensions
n = len(matrix)

# Create rotated matrix
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Fill rotated matrix
for i in range(n):
    for j in range(n):
        rotated[j][n-1-i] = matrix[i][j]

# Format output string
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"

print(output.strip())