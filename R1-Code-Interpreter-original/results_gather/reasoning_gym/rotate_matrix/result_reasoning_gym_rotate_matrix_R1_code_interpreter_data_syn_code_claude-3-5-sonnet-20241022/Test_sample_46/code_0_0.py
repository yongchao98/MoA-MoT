# Create the input matrix
matrix = [
    [1, 4, 4, 6, 0],
    [2, 6, 3, 7, 6],
    [8, 0, 3, 2, 5],
    [6, 5, 6, 2, 2],
    [5, 5, 0, 3, 1]
]

# For 270° clockwise (or 90° counterclockwise)
n = len(matrix)
rotated = [[0 for _ in range(n)] for _ in range(n)]

for i in range(n):
    for j in range(n):
        rotated[n-1-j][i] = matrix[i][j]

# Format the output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
output = output.strip()

print("<<<" + output + ">>>")