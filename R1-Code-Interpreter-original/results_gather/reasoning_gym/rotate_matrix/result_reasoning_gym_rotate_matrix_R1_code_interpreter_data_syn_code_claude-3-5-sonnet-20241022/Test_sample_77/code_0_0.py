# Input matrix
matrix = [
    [3, 7, 7, 8, 0, 1],
    [3, 8, 7, 4, 9, 8],
    [1, 3, 7, 8, 4, 6],
    [2, 4, 8, 8, 5, 1],
    [8, 2, 0, 5, 9, 3],
    [2, 9, 2, 6, 7, 2]
]

n = len(matrix)
# Create rotated matrix
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Rotate 180 degrees
for i in range(n):
    for j in range(n):
        rotated[i][j] = matrix[n-1-i][n-1-j]

# Format output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
output = output.strip()

print("<<<" + output + ">>>")