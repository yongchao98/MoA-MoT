# Input matrix
matrix = [[5, 9],
         [9, 0]]

# For 270 degrees clockwise (or 90 degrees counterclockwise)
# We can do this in one step rather than three 90-degree rotations
n = len(matrix)
rotated = [[0 for _ in range(n)] for _ in range(n)]

for i in range(n):
    for j in range(n):
        rotated[n-1-j][i] = matrix[i][j]

# Format output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
output = output.strip()

print(output)