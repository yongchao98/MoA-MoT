# Input matrix as string
matrix_str = """3 3 0 3 7
6 5 4 4 7
9 6 4 9 1
0 7 3 7 1
7 8 6 5 0"""

# Convert string to matrix
matrix = [list(map(int, row.split())) for row in matrix_str.strip().split('\n')]

# Rotate matrix 90 degrees clockwise
n = len(matrix)
rotated = [[0 for _ in range(n)] for _ in range(n)]

for i in range(n):
    for j in range(n):
        rotated[j][n-1-i] = matrix[i][j]

# Format output
result = []
for row in rotated:
    result.append(' '.join(map(str, row)))
result = '\n'.join(result)

print(f"<<<{result}>>>")