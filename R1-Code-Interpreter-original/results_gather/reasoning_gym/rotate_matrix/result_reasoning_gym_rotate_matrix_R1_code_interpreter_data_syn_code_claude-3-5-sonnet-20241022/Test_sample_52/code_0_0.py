# Input matrix
matrix = '''0 7 0 8 5 3
7 7 8 5 9 2
8 4 8 7 6 5
9 3 2 7 0 8
0 2 0 7 8 4
8 1 0 3 3 6'''

# Convert string to 2D list
matrix = [list(map(int, row.split())) for row in matrix.split('\n')]
n = len(matrix)

# Create rotated matrix
rotated = [[0 for _ in range(n)] for _ in range(n)]

# Rotate 180 degrees
for i in range(n):
    for j in range(n):
        rotated[i][j] = matrix[n-1-i][n-1-j]

# Format output
output = '\n'.join(' '.join(str(x) for x in row) for row in rotated)
print(f"<<<{output}>>>")